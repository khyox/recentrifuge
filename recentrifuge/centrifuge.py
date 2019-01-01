"""
Functions directly related with Centrifuge (or Kraken).

"""

import collections as col
import io
import os
from math import log10
from statistics import mean
from typing import Tuple, Counter, Dict, List, Set

from Bio import SeqIO

from recentrifuge.config import Filename, Id, Score, Scoring
from recentrifuge.config import gray, red, green, yellow
from recentrifuge.rank import Rank
from recentrifuge.stats import SampleStats

# Centrifuge specific constants
UNCLASSIFIED: Id = Id('0')


def read_report(report_file: str) -> Tuple[str, Counter[Id],
                                           Dict[Id, Rank]]:
    """
    Read Centrifuge/Kraken report file

    Args:
        report_file: report file name

    Returns:
        log string, abundances counter, taxlevel dict

    """
    # TODO: Discontinued method, to be erased in a future release
    output: io.StringIO = io.StringIO(newline='')
    abundances: Counter[Id] = col.Counter()
    level_dic = {}
    output.write(f'\033[90mLoading report file {report_file}...\033[0m')
    try:
        with open(report_file, 'r') as file:
            for report_line in file:
                _, _, taxnum, taxlev, _tid, _ = report_line.split('\t')
                tid = Id(_tid)
                abundances[tid] = int(taxnum)
                level_dic[tid] = Rank.centrifuge(taxlev)
    except KeyboardInterrupt:
        print(gray(' User'), yellow('interrupted!'))
        raise
    except Exception:
        print(red('ERROR!'), 'Cannot read "' + report_file + '"')
        raise
    else:
        output.write('\033[92m OK! \033[0m\n')
    return output.getvalue(), abundances, level_dic


def read_output(output_file: Filename,
                scoring: Scoring = Scoring.SHEL,
                minscore: Score = None,
                ) -> Tuple[str, SampleStats,
                           Counter[Id], Dict[Id, Score]]:
    """
    Read Centrifuge output file

    Args:
        output_file: output file name
        scoring: type of scoring to be applied (see Scoring class)
        minscore: minimum confidence level for the classification

    Returns:
        log string, statistics, abundances counter, scores dict

    """
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
    try:
        with open(output_file, 'r') as file:
            file.readline()  # discard header
            for output_line in file:
                try:
                    _, _, _tid, _score, _, _, _length, *_ = output_line.split(
                        '\t')
                except ValueError:
                    print(yellow('Failure'), 'parsing line elements:'
                                             f' {output_line} in {output_file}'
                                             '. Ignoring line!')
                    last_error_read = num_read + 1
                    num_errors += 1
                    continue
                tid = Id(_tid)
                try:
                    # From Centrifuge score get "single hit equivalent length"
                    shel = Score(float(_score) ** 0.5 + 15)
                    length = int(_length)
                except ValueError:
                    print(yellow('Failure'), f'parsing score ({_score}) for ',
                          f'query length {_length} for taxid {_tid}',
                          f'in {output_file}. Ignoring line!')
                    last_error_read = num_read + 1
                    num_errors += 1
                    continue
                num_read += 1
                nt_read += length
                if tid == UNCLASSIFIED:  # Just count unclassified reads
                    num_uncl += 1
                    continue
                else:
                    taxids.add(tid)  # Save all the tids of classified reads
                if minscore is not None and shel < minscore:
                    continue  # Ignore read if low confidence
                try:
                    all_scores[tid].append(shel)
                except KeyError:
                    all_scores[tid] = [shel, ]
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
        minscore=minscore, nt_read=nt_read, scores=all_scores, lens=all_length,
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
    output.write(gray('  Scores: min = ') + f'{stat.sco.mini:.1f}' +
                 gray(', max = ') + f'{stat.sco.maxi:.1f}' +
                 gray(', avr = ') + f'{stat.sco.mean:.1f}\n')
    output.write(gray('  Length: min = ') + f'{stat.len.mini}' +
                 gray(', max = ') + f'{stat.len.maxi}' +
                 gray(', avr = ') + f'{stat.len.mean}\n')
    output.write(gray('  TaxIds: by classifier = ') + f'{stat.tid.clas}'
                 + gray(', by filter = ') + f'{stat.tid.filt}\n')
    # Select score output
    out_scores: Dict[Id, Score]
    if scoring is Scoring.SHEL:
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
        print(red('ERROR!'), f' Centrifuge: Unsupported Scoring "{scoring}"')
        raise Exception('Unsupported scoring')
    # Return
    return output.getvalue(), stat, counts, out_scores


def analize_outputs(outputs: List[Filename]) -> None:
    """
    Check several Centrifuge output files.

    Args:
        outputs: is a list of Centrifuge output files.
    """
    for cfgoutfile in outputs:
        with open(cfgoutfile, 'rU') as lmo_handle:
            print(f'For cfg output {cfgoutfile}...')
            seqs_lst = list(SeqIO.parse(lmo_handle, 'centrifuge'))
            print(f'- There are {len(seqs_lst)} seqs processed.')
            scores = [seq.annotations['score'] for seq in seqs_lst if
                      seq.annotations['score'] > 0]
            print(f'- Scores: max={max(scores)}, min={min(scores)}')


def select_centrifuge_inputs(outputs: List[Filename],
                             ext: str = '.out') -> None:
    """Centrifuge output files processing specific stuff"""
    dir_name = outputs[0]
    outputs.clear()
    with os.scandir(dir_name) as dir_entry:
        for fil in dir_entry:
            if not fil.name.startswith('.') and fil.name.endswith(ext):
                if dir_name != '.':
                    outputs.append(Filename(os.path.join(dir_name, fil.name)))
                else:  # Avoid sample names starting with just the dot
                    outputs.append(Filename(fil.name))
    outputs.sort()
    print(gray(f'Centrifuge {ext} files to analyze:'), outputs)
