"""
Functions directly related with Kraken.

"""

import collections as col
import io
import os
import re
from math import log10
from statistics import mean
from typing import Tuple, Counter, Dict, List, Set, Any, Union, TextIO, IO

from recentrifuge.config import Filename, Id, Score, Scoring, GO_ROOT
from recentrifuge.config import gray, red, green, yellow, blue, magenta
from recentrifuge.stats import SampleStats

# KRAKEN specific constants
UNCLASSIFIED: str = 'U'
K_MER_SIZE: int = 35  # Default k-mer size for Kraken


def open_compressed_and_uncompressed(filename: Filename
                                     ) -> Union[TextIO, IO[Any]]:
    """Aux method to deal with kraken compressed output"""
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        import gzip
        return gzip.open(filename, mode='rt')
    elif ext == '.bz2':
        import bz2
        return bz2.open(filename, mode='rt')
    else:
        return open(filename, mode='rt')


def read_kraken_output(output_file: Filename,
                      scoring: Scoring = Scoring.KRAKEN,
                      minscore: Score = None,
                      ) -> Tuple[str, SampleStats,
                           Counter[Id], Dict[Id, Score]]:
    """
    Read Kraken output file

    Args:
        output_file: output file name
        scoring: type of scoring to be applied (see Scoring class)
        minscore: minimum confidence level for the classification

    Returns:
        log string, statistics, abundances counter, scores dict

    """
    output: io.StringIO = io.StringIO(newline='')
    all_scores: Dict[Id, List[Score]] = {}
    all_kmerel: Dict[Id, List[Score]] = {}
    all_length: Dict[Id, List[int]] = {}
    taxids: Set[Id] = set()
    num_read: int = 0
    nt_read: int = 0
    num_uncl: int = 0
    last_error_read: int = -1  # Number of read of the last error
    num_errors: int = 0  # Number or reads discarded due to error
    output.write(gray(f'Loading output file {output_file}... '))
    try:
        with open_compressed_and_uncompressed(output_file) as file:
            check_format: bool = True
            are_proteins: bool = False
            for raw_line in file:
                if check_format:  # Check number of cols in first line
                    header = raw_line.split('\t')
                    if len(header) != 5:
                        print(red('\nERROR! ') + 'Kraken output format of ',
                              yellow(f'"{output_file}"'), 'not supported.')
                        print(magenta('Expected:'),
                              'C/U, seqID, go/taxID, length, list of mappings')
                        print(magenta('Found:'), '\t'.join(header), end='')
                        print(blue('HINT:'),
                              'Use Kraken or Kraken2 direct output.')
                        raise Exception('Unsupported file format. Aborting.')
                try:
                    output_line = raw_line.strip()
                    (_clas, _label, _myid, _length,
                     _maps) = output_line.split('\t')
                except ValueError:
                    print(yellow('Failure'), 'parsing line elements:'
                                             f' {output_line} in {output_file}'
                                             '. Ignoring line!')
                    last_error_read = num_read + 1
                    num_errors += 1
                    continue
                try:
                    maps: Union[List[str], Set[str]] = _maps.split()
                    if check_format:
                        if '-:-' in maps:  # Results from protein classification
                            are_proteins = True
                        check_format = False  # Finished checking format
                    length: int = sum(map(int, _length.split('|')))
                    num_read += 1
                    nt_read += length
                    if _clas == UNCLASSIFIED:  # Just count unclassified reads
                        num_uncl += 1
                        continue
                    myid: Id
                    if are_proteins:  # TODO: Check again when krk2 DB is right!
                        myid = Id('GO:' + _myid.zfill(7))
                        if myid == 'GO:0000001':  # TODO: Change virtual root of GO in DB build
                            myid = GO_ROOT
                    else:
                        try:
                            myid = Id(str(int(_myid)))
                        except ValueError:
                            myid = Id(re.search('\(taxid (\d+)\)', _myid).group(1))
                    if are_proteins:
                        # TODO: Check for better approaches (6 frames problem)
                        maps = set(maps)  # Remove duplicated frame entries
                        maps.remove('-:-')
                    else:  # Results from taxonomic classification
                        try:  # try to remove paired read separator
                            maps.remove('|:|')
                        except ValueError:
                            pass
                    mappings: Counter[Id] = col.Counter()
                    for pair in maps:
                        couple: List[str] = pair.split(':')
                        if are_proteins:
                            mappings[Id(
                                'GO:' + couple[0].zfill(7))] += int(couple[1])
                        else:
                            mappings[Id(couple[0])] += int(couple[1])
                    # From Kraken score get "single hit equivalent length"
                    shel: Score = Score(mappings[myid] + K_MER_SIZE)
                    score: Score = Score(mappings[myid] / sum(mappings.values())
                                         * 100)  # % relative to all k-mers
                except (ValueError, IndexError):
                    print(yellow('Failure'), 'parsing line elements:'
                                             f' {output_line} in {output_file}'
                                             '. Ignoring line!')
                    last_error_read = num_read + 1
                    num_errors += 1
                    continue
                else:
                    taxids.add(myid)  # Save all the tids of classified reads
                if minscore is not None:  # Decide if ignore read if low score
                    if scoring is Scoring.KRAKEN:
                        if score < minscore:
                            continue
                    else:
                        if shel < minscore:
                            continue
                try:
                    all_scores[myid].append(shel)
                except KeyError:
                    all_scores[myid] = [shel, ]
                try:
                    all_kmerel[myid].append(score)
                except KeyError:
                    all_kmerel[myid] = [score, ]
                try:
                    all_length[myid].append(length)
                except KeyError:
                    all_length[myid] = [length, ]
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
        minscore=minscore, nt_read=nt_read, lens=all_length,
        scores=all_scores, scores2=all_kmerel,
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
    output.write(gray('  Scores SHEL: min = ') + f'{stat.sco.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco.mean:.1f}\n')
    output.write(gray('  Coverage(%): min = ') + f'{stat.sco2.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco2.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco2.mean:.1f}\n')
    output.write(gray('  Read length: min = ') + f'{stat.len.mini},' +
                 gray(' max = ') + f'{stat.len.maxi},' +
                 gray(' avr = ') + f'{stat.len.mean}\n')
    output.write(gray('  TaxIds: by classifier = ') + f'{stat.tid.clas}'
                 + gray(', by filter = ') + f'{stat.tid.filt}\n')
    # Select score output
    out_scores: Dict[Id, Score]
    if scoring is Scoring.SHEL:
        out_scores = {tid: Score(mean(all_scores[tid])) for tid in all_scores}
    elif scoring is Scoring.KRAKEN:
        out_scores = {tid: Score(mean(all_kmerel[tid])) for tid in all_kmerel}
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
        print(red('ERROR!'), f'kraken: Unsupported Scoring "{scoring}"')
        raise Exception('Unsupported scoring')
    # Return
    return output.getvalue(), stat, counts, out_scores


def select_kraken_inputs(krakens: List[Filename],
                         ext: str = '.krk') -> None:
    """Search for Kraken files to analyze"""
    dir_name = krakens[0]
    krakens.clear()
    with os.scandir(dir_name) as dir_entry:
        for fil in dir_entry:
            if not fil.name.startswith('.') and fil.name.endswith(ext):
                if dir_name != '.':
                    krakens.append(Filename(os.path.join(dir_name, fil.name)))
                else:  # Avoid sample names starting with just the dot
                    krakens.append(Filename(fil.name))
    krakens.sort()
    print(gray(f'Kraken {ext} files to analyze:'), krakens)
