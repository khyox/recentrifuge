"""
Functions directly related with a generic taxonomic classifier.

"""

import collections as col
import io
from enum import Enum
from math import log10
from statistics import mean
from typing import Tuple, Counter, Dict, List, Set

from recentrifuge.config import Filename, Id, Score, Scoring
from recentrifuge.config import gray, red, green, yellow, blue
from recentrifuge.stats import SampleStats


class GenericType(Enum):
    """Enumeration with options for the type of supported generic file."""
    CSV = 0  # Comma Separated Values
    TSV = 1  # Tab Separated Values
    SSV = 2  # Space Separated Values

    def __str__(self):
        return f'{str(self.name)}'


class GenericFormat(object):
    """Class containing the format of a a generic classifier file."""

    MIN_COLS: int = 3  # Minimum number of columns: TaxID, LENgth, SCOre
    MIN_FIELDS: int = MIN_COLS + 2  # Plus TYPe and UNClassified

    def __init__(self, frmt: str):

        def print_error(specifier):
            """GenericFormat constructor: print an informative error message"""
            print(red('ERROR!'), 'Generic --format string malformed:',
                  blue(specifier), '\n\tPlease rerun with --help for details.')

        blocks: List[str] = [fld.strip() for fld in frmt.split(',')]
        if len(blocks) < self.MIN_FIELDS:
            print_error(f'Wrong number of fields (expected {self.MIN_FIELDS} '
                        f'found {len(blocks)}).')
            exit(2)
        try:
            fmt: Dict[str, str] = {pair.split(':')[0].strip():
                                   pair.split(':')[1].strip()
                                   for pair in blocks}
        except IndexError:
            print_error('All fields need ":" separator.')
            exit(2)
        # Populate fields
        try:
            typ = fmt['TYP']
        except KeyError:
            print_error('TYPe field is mandatory.')
            exit(2)
        try:
            self.typ: GenericType = GenericType[typ.upper()]
        except KeyError:
            print_error('Unknown file TYPe, valid options are ' +
                        ' or '.join([str(t) for t in GenericType]))
            exit(2)
        try:
            self.tid: int = int(fmt['TID'])
        except KeyError:
            print_error('TaxID field is mandatory.')
            exit(2)
        except ValueError:
            print_error('TaxID field is an integer number of column.')
            exit(2)
        try:
            self.len: int = int(fmt['LEN'])
        except KeyError:
            print_error('LENgth field is mandatory.')
            exit(2)
        except ValueError:
            print_error('LENgth field is an integer number of column.')
            exit(2)
        try:
            self.sco: int = int(fmt['SCO'])
        except KeyError:
            print_error('SCOre field is mandatory.')
            exit(2)
        except ValueError:
            print_error('SCOre field is an integer number of column.')
            exit(2)
        try:
            self.unc: Id = Id(fmt['UNC'])
        except KeyError:
            print_error('UNClassified field is mandatory.')
            exit(2)
        # Check columns are different
        if (self.tid == self.len or self.tid == self.sco
                or self.len == self.sco):
            print_error('Different fields need different columns.')
            exit(2)
        # Check column numbers are positive
        if self.tid < 1 or self.len < 1 or self.sco < 1:
            print_error('Columns numbers should be positive integers.')
            exit(2)

    def __str__(self):
        return (f'Generic format = TYP:{self.typ}, TID:{self.tid}, '
                f'LEN:{self.len}, SCO:{self.sco}, UNC:{self.unc}.')


def read_generic_output(output_file: Filename,
                        scoring: Scoring = Scoring.GENERIC,
                        minscore: Score = None,
                        genfmt: GenericFormat = None
                        ) -> Tuple[str, SampleStats,
                                   Counter[Id], Dict[Id, Score]]:
    """
    Read an output file from a generic classifier

    Args:
        output_file: output file name
        scoring: type of scoring to be applied (see Scoring class)
        minscore: minimum confidence level for the classification
        genfmt: GenericFormat object specifying the files format

    Returns:
        log string, statistics, abundances counter, scores dict

    """
    # Initialization of variables
    output: io.StringIO = io.StringIO(newline='')
    all_scores: Dict[Id, List[Score]] = {}
    all_length: Dict[Id, List[int]] = {}
    taxids: Set[Id] = set()
    num_read: int = 0
    nt_read: int = 0
    num_uncl: int = 0
    last_error_read: int = -1  # Number of read of the last error
    num_errors: int = 0  # Number or reads discarded due to error
    output.write(gray(f'Loading output file {output_file}... '))
    # Check format
    if not isinstance(genfmt, GenericFormat):
        raise Exception(red('\nERROR!'),
                        'Missing GenericFormat when reading a generic output.')
    try:
        with open(output_file, 'r') as file:
            # Main loop processing each file line
            for raw_line in file:
                raw_line = raw_line.strip(' \n\t')
                splitting: str
                if genfmt.typ is GenericType.CSV:
                    splitting = ','
                elif genfmt.typ is GenericType.TSV:
                    splitting = '\t'
                elif genfmt.typ is GenericType.SSV:
                    splitting = ' '
                else:
                    raise Exception(f'ERROR! Unknown GenericType {genfmt.typ}')
                output_line: List[str] = raw_line.split(splitting)
                if len(output_line) < GenericFormat.MIN_COLS:
                    if num_read == 0 and last_error_read < 0:
                        last_error_read = 0
                        print(yellow('Warning!'), 'Skipping header of '
                                                  f'{output_file}')
                        continue  # Not account for the header as an error
                    raise Exception(
                        red('\nERROR!') + ' Line ' + yellow(f'{output_line}')
                        + '\n\tin ' + yellow(f'{output_file}') + ' has < '
                        + blue(f'{GenericFormat.MIN_COLS}') + ' required '
                        + 'columns.\n\tPlease check the file.')
                try:
                    tid: Id = Id(output_line[genfmt.tid-1].strip(' "'))
                    length: int = int(output_line[genfmt.len-1].strip(' "'))
                    if tid == genfmt.unc:  # Avoid read score for unclass reads
                        num_read += 1
                        nt_read += length
                        num_uncl += 1
                        continue
                    score: Score = Score(
                        float(output_line[genfmt.sco-1].strip(' "')))
                except ValueError:
                    if num_read == 0 and last_error_read < 0:
                        last_error_read = 0
                        print(yellow('Warning!'), 'Skipping header of '
                                                  f'{output_file}')
                        continue  # Not account for the header as a failure
                    print(yellow('Failure'), 'parsing line elements:'
                                             f' {output_line} in {output_file}'
                                             '. Ignoring line!')
                    last_error_read = num_read + 1
                    num_errors += 1
                    if num_read > 100 and num_errors > 0.5 * num_read:
                        print(red('ERROR!'),
                              'Unreliable file processing: rate of problematic'
                              f' reads is {num_errors/num_read*100:_d}, beyond'
                              ' 50%, after 100 reads. Please check the format '
                              f'of the file "{output_file}".')
                        raise
                    else:
                        continue
                num_read += 1
                nt_read += length
                taxids.add(tid)  # Save all the tids of classified reads
                if minscore is not None and score < minscore:
                    continue  # Discard read if low confidence
                try:
                    all_scores[tid].append(score)
                except KeyError:
                    all_scores[tid] = [score, ]
                try:
                    all_length[tid].append(length)
                except KeyError:
                    all_length[tid] = [length, ]
    except FileNotFoundError:
        raise Exception(red('\nERROR! ') + f'Cannot read "{output_file}"')
    if last_error_read == num_read + 1:  # Check error in last line: truncated!
        print(yellow('Warning!'), f'{output_file} seems truncated!')
    counts: Counter[Id] = col.Counter({tid: len(all_scores[tid])
                                       for tid in all_scores})
    output.write(green('OK!\n'))
    if num_read == 0:
        raise Exception(red('\nERROR! ')
                        + f'Cannot read any sequence from "{output_file}"')
    filt_seqs: int = sum([len(scores) for scores in all_scores.values()])
    if filt_seqs == 0:
        raise Exception(red('\nERROR! ') + 'No sequence passed the filter!')
    # Get statistics
    stat: SampleStats = SampleStats(
        minscore=minscore, nt_read=nt_read, lens=all_length, scores=all_scores,
        seq_read=num_read, seq_unclas=num_uncl, seq_filt=filt_seqs,
        tid_clas=len(taxids)
    )
    # Output statistics
    if num_errors:
        output.write(gray('  Seqs fail: ') + red(f'{num_errors:_d}\t') +
                     gray('(Last error in read ') + red(f'{last_error_read}') +
                     gray(')\n'))
    output.write(gray('  Seqs read: ') + f'{stat.seq.read:_d}\t' + gray('[')
                 + f'{stat.nt_read}' + gray(']\n'))
    output.write(gray('  Seqs clas: ') + f'{stat.seq.clas:_d}\t' + gray('(') +
                 f'{stat.get_unclas_ratio():.2%}' + gray(' unclassified)\n'))
    output.write(gray('  Seqs pass: ') + f'{stat.seq.filt:_d}\t' + gray('(') +
                 f'{stat.get_reject_ratio():.2%}' + gray(' rejected)\n'))
    output.write(gray('  Scores: min = ') + f'{stat.sco.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco.mean:.1f}\n')
    output.write(gray('  Read length: min = ') + f'{stat.len.mini},' +
                 gray(' max = ') + f'{stat.len.maxi},' +
                 gray(' avr = ') + f'{stat.len.mean}\n')
    output.write(gray('  TaxIds: by classifier = ') + f'{stat.tid.clas}'
                 + gray(', by filter = ') + f'{stat.tid.filt}\n')
    # Select score output
    out_scores: Dict[Id, Score]
    if scoring is Scoring.GENERIC:
        out_scores = {tid: Score(mean(all_scores[tid])) for tid in all_scores}
    elif scoring is Scoring.LENGTH:
        out_scores = {tid: Score(mean(all_length[tid])) for tid in all_length}
    elif scoring is Scoring.LOGLENGTH:
        out_scores = {tid: Score(log10(mean(all_length[tid])))
                      for tid in all_length}
    elif scoring is Scoring.NORMA:
        scores: Dict[Id, Score] = {tid: Score(mean(all_scores[tid]))
                                   for tid in all_scores}
        lengths: Dict[Id, Score] = {tid: Score(mean(all_length[tid]))
                                    for tid in all_length}
        out_scores = {tid: Score(scores[tid] / lengths[tid] * 100)
                      for tid in scores}
    else:
        raise Exception(red('\nERROR!'),
                        f'Generic: Unsupported Scoring "{scoring}"')
    # Return
    return output.getvalue(), stat, counts, out_scores
