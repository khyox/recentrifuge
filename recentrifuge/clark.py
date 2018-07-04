"""
Functions directly related with CLARK_C(S).

"""

import io
import os
from math import log10
from statistics import mean
from typing import Tuple, Counter, Dict, List

from recentrifuge.config import Filename, Id, Score, Scoring
from recentrifuge.config import gray, red, green, yellow
from recentrifuge.stats import SampleStats

# CLARK_C specific constants
UNCLASSIFIED: Id = Id('NA')


def read_clark_output(output_file: Filename,
                      scoring: Scoring = Scoring.CLARK_C,
                      minscore: Score = None,
                      ) -> Tuple[str, SampleStats,
                           Counter[Id], Dict[Id, Score]]:
    """
    Read CLARK_C(S) output file

    Args:
        output_file: output file name
        scoring: type of scoring to be applied (see Scoring class)
        minscore: minimum confidence level for the classification

    Returns:
        log string, statistics, abundances counter, scores dict

    """
    output: io.StringIO = io.StringIO(newline='')
    all_scores: Dict[Id, List[Score]] = {}
    all_confs: Dict[Id, List[Score]] = {}
    all_gammas: Dict[Id, List[Score]] = {}
    all_length: Dict[Id, List[int]] = {}
    num_read: int = 0
    nt_read: int = 0
    num_uncl: int = 0
    error_read: int = -1
    output.write(gray(f'Loading output file {output_file}... '))
    try:
        with open(output_file, 'r') as file:
            # Check number of cols in header
            header = file.readline().split(',')
            if len(header) != 8:
                print(red('\nERROR! ') +
                      f'CLARK_C output format of "{output_file}" not supported.')
                print('Expected: ID,Length,Gamma,1st,score1,2nd,score2,conf')
                raise Exception('Unsupported file format. Aborting.')
            for raw_line in file:
                try:
                    output_line = raw_line.strip()
                    (_label, _length, _gamma, _tid1, _score1, _tid2, _score2,
                     _conf) = output_line.split(',')
                except ValueError:
                    print(yellow('Error'), f' parsing line: ({output_line}) '
                                           f'in {output_file}. Ignoring line!')
                    error_read = num_read + 1
                    continue
                try:
                    length: int = int(_length)
                    gamma: Score = Score(float(_gamma))
                    tid1: Id = Id(_tid1)
                    score1: Score = Score(float(_score1))
                    tid2: Id = Id(_tid2)
                    score2: Score = Score(float(_score2))
                    conf: Score = Score(float(_conf))
                except ValueError:
                    print(yellow('Error'), f' parsing elements of'
                                           f' line: ({output_line}) '
                                           f'in {output_file}. Ignoring line!')
                    error_read = num_read + 1
                    continue
                num_read += 1
                nt_read += length
                # Select tid and score between CLARK_C assignments 1 and 2
                tid: Id = tid1
                score: Score = score1
                if tid1 == UNCLASSIFIED:
                    if tid2 == UNCLASSIFIED:  # Just count unclassified reads
                        num_uncl += 1
                        continue
                    else:  # Majority of read unclassified
                        tid = tid2
                        score = score2
                        conf = Score(1 - conf)  # Get CLARK_C's h2/(h1+h2)
                # From CLARK_C(S) score get "single hit equivalent length"
                shel: Score = Score(score)
                if minscore is not None:  # Decide if ignore read if low score
                    if scoring is Scoring.CLARK_C:
                        if conf < minscore:
                            continue
                    elif scoring is Scoring.CLARK_G:
                        if gamma < minscore:
                            continue
                    else:
                        if shel < minscore:
                            continue
                try:
                    all_scores[tid].append(shel)
                except KeyError:
                    all_scores[tid] = [shel, ]
                try:
                    all_confs[tid].append(conf)
                except KeyError:
                    all_confs[tid] = [conf, ]
                try:
                    all_gammas[tid].append(gamma)
                except KeyError:
                    all_gammas[tid] = [gamma, ]
                try:
                    all_length[tid].append(length)
                except KeyError:
                    all_length[tid] = [length, ]

    except FileNotFoundError:
        raise Exception(red('\nERROR! ') + f'Cannot read "{output_file}"')
    if error_read == num_read + 1:  # Check if error in last line: truncated!
        print(yellow('Warning!'), f'{output_file} seems truncated!')
    counts: Counter[Id] = Counter({tid: len(all_scores[tid])
                                   for tid in all_scores})
    output.write(green('OK!\n'))
    if num_read == 0:
        raise Exception(red('\nERROR! ')
                        + f'Cannot read any sequence from"{output_file}"')
    filt_seqs: int = sum([len(scores) for scores in all_scores.values()])
    if filt_seqs == 0:
        raise Exception(red('\nERROR! ') + 'No sequence passed the filter!')
    # Get statistics
    stat: SampleStats = SampleStats(
        minscore=minscore, nt_read=nt_read, lens=all_length,
        scores=all_scores, scores2=all_confs, scores3=all_gammas,
        seq_read=num_read, seq_unclas=num_uncl, seq_filt=filt_seqs
    )
    # Output statistics
    output.write(gray('  Seqs read: ') + f'{stat.seq.read:_d}\t' + gray('[')
                 + f'{stat.nt_read}' + gray(']\n'))
    output.write(gray('  Seqs clas: ') + f'{stat.seq.clas:_d}\t' + gray('(') +
                 f'{stat.get_unclas_ratio():.2%}' + gray(' unclassified)\n'))
    output.write(gray('  Seqs pass: ') + f'{stat.seq.filt:_d}\t' + gray('(') +
                 f'{stat.get_reject_ratio():.2%}' + gray(' rejected)\n'))
    output.write(gray('  Hit (score): min = ') + f'{stat.sco.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco.mean:.1f}\n')
    output.write(gray('  Conf. score: min = ') + f'{stat.sco2.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco2.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco2.mean:.1f}\n')
    output.write(gray('  Gamma score: min = ') + f'{stat.sco3.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco3.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco3.mean:.1f}\n')
    output.write(gray('  Read length: min = ') + f'{stat.len.mini},' +
                 gray(' max = ') + f'{stat.len.maxi},' +
                 gray(' avr = ') + f'{stat.len.mean}\n')
    output.write(f'  {stat.num_taxa}' + gray(f' taxa with assigned reads\n'))
    # Select score output
    out_scores: Dict[Id, Score]
    if scoring is Scoring.SHEL:
        out_scores = {tid: Score(mean(all_scores[tid])) for tid in all_scores}
    elif scoring is Scoring.CLARK_C:
        out_scores = {tid: Score(mean(all_confs[tid])) for tid in all_confs}
    elif scoring is Scoring.CLARK_G:
        out_scores = {tid: Score(mean(all_gammas[tid])) for tid in all_gammas}
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
        print(red('ERROR!'), f'clark: Unsupported Scoring "{scoring}"')
        raise Exception('Unsupported scoring')
    # Return
    return output.getvalue(), stat, counts, out_scores


def select_clark_inputs(clarks: List[Filename],
                        ext: str = '.csv') -> None:
    """CLARK_C output files processing specific stuff"""
    dir_name = clarks[0]
    clarks.clear()
    with os.scandir(dir_name) as dir_entry:
        for fil in dir_entry:
            if not fil.name.startswith('.') and fil.name.endswith(ext):
                if dir_name != '.':
                    clarks.append(Filename(os.path.join(dir_name, fil.name)))
                else:  # Avoid sample names starting with just the dot
                    clarks.append(Filename(fil.name))
    clarks.sort()
    print(gray(f'CLARK_C {ext} files to analyze:'), clarks)