"""
Functions directly related with Gene Ontology.

"""

import collections as col
import io
import os
import sys
import time
from math import log10, isclose
from statistics import mean
from typing import Tuple, Counter, Set, Dict, List, Optional, Union

from recentrifuge.config import Filename, Id, Score, Scoring, Sample
from recentrifuge.config import Err, NO_SCORE
from recentrifuge.config import gray, red, green, yellow, blue
from recentrifuge.kraken import read_kraken_output
from recentrifuge.ontology import Ontology
from recentrifuge.stats import SampleStats
from recentrifuge.trees import TaxTree, SampleDataById
from recentrifuge.rank import Ranks


def read_blast_go(output_file: Filename,
                  scoring: Scoring = Scoring.SHEL,
                  minscore: Score = None,
                  ) -> Tuple[str, SampleStats,
                       Counter[Id], Dict[Id, Score]]:
    """
    Read GO output file

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
    num_read: int = 0
    nt_read: int = 0
    num_uncl: int = 0  # TODO: In the future, it could account for GI->GO fails
    error_read: int = None
    output.write(gray(f'Loading GO data file {output_file}... '))
    try:
        with open(output_file, 'r') as file:
            file.readline()  # discard header
            for output_line in file:
                try:
                    _tid, _gi, _go, _score, _evalue, *_ = output_line.rstrip(
                        ).split('\t')
                except ValueError:
                    print(red('Error'), f'parsing line: ({output_line}) '
                                        f'in {output_file}. Ignoring line!')
                    error_read = num_read + 1
                    continue
                # tid = Id(_tid)
                gid: Id = Id(_go)
                try:
                    # From Centrifuge score get "single hit equivalent length"
                    evalue: float = float(_evalue)
                    score: Score
                    if isclose(evalue, 1e-308, abs_tol=1e-308):
                        score = Score(308)
                    else:
                        score = Score(-log10(evalue))
                    length = int(_score)
                except ValueError:
                    print(red('Error'), f'parsing eval ({_evalue}) for query',
                          f'score ({_score}) for GO {_go}',
                          f'in {output_file}. Ignoring line!')
                    continue
                num_read += 1
                nt_read += length
                if minscore is not None and score < minscore:
                    continue  # Ignore read if low confidence
                try:
                    all_scores[gid].append(score)
                except KeyError:
                    all_scores[gid] = [score, ]
                try:
                    all_length[gid].append(length)
                except KeyError:
                    all_length[gid] = [length, ]
    except FileNotFoundError:
        raise Exception(red('\nERROR! ') + f'Cannot read "{output_file}"')
    if error_read == num_read + 1:  # Check if error in last line: truncated!
        print(yellow('Warning!'), f'{output_file} seems truncated!')
    counts: Counter[Id] = col.Counter({tid: len(all_scores[tid])
                                   for tid in all_scores})
    output.write(green('OK!\n'))
    if num_read == 0:
        raise Exception(red('\nERROR! ')
                        + f'Cannot read any sequence from"{output_file}"')
    filt_seqs: int = sum([len(scores) for scores in all_scores.values()])
    if filt_seqs == 0:
        raise Exception(red('\nERROR! ') + f' In "{output_file}" '
                                           'no sequence passed the filter!')
    # Get statistics
    stat: SampleStats = SampleStats(
        minscore=minscore, nt_read=nt_read, scores=all_scores, lens=all_length,
        seq_read=num_read, seq_unclas=num_uncl, seq_filt=filt_seqs
    )
    # Output statistics
    output.write(gray('  Annotations read: ') + f'{stat.seq.read:_d}\t'
                 + gray('[') + f'{stat.nt_read}' + gray(']\n'))
# TODO: In the future, it could account for GI->GO fails
#    output.write(gray('  Seqs clas: ') + f'{stat.seq.clas:_d}\t' + gray('(') +
#                  f'{stat.get_unclas_ratio():.2%}' + gray(' unclassified)\n'))
    output.write(gray('  Annotations pass: ') + f'{stat.seq.filt:_d}\t' +
                 gray('(') +f'{stat.get_reject_ratio():.2%}' +
                 gray(' rejected)\n'))
    output.write(gray('  Scores: min = ') + f'{stat.sco.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco.mean:.1f}\n')
    output.write(gray('  Length: min = ') + f'{stat.len.mini},' +
                 gray(' max = ') + f'{stat.len.maxi},' +
                 gray(' avr = ') + f'{stat.len.mean}\n')
    output.write(gray('  GOids: by classifier = ') + f'{stat.tid.clas}'
                 + gray(', by filter = ') + f'{stat.tid.filt}\n')

    # Select score output
    out_scores: Dict[Id, Score]
    if scoring is Scoring.SHEL:
        out_scores = {gid: Score(mean(all_scores[gid])) for gid in all_scores}
    elif scoring is Scoring.LENGTH:
        out_scores = {gid: Score(mean(all_length[gid])) for gid in all_length}
    elif scoring is Scoring.LOGLENGTH:
        out_scores = {gid: Score(log10(mean(all_length[gid])))
                      for gid in all_length}
    elif scoring is Scoring.NORMA:
        scores: Dict[Id, Score] = {gid: Score(mean(all_scores[gid]))
                                   for gid in all_scores}
        lengths: Dict[Id, Score] = {gid: Score(mean(all_length[gid]))
                                    for gid in all_length}
        out_scores = {gid: Score(scores[gid] / lengths[gid] * 100)
                      for gid in scores}
    else:
        raise Exception(f'\n\033[91mERROR!\033[0m Unknown Scoring "{scoring}"')
    # Return
    return output.getvalue(), stat, counts, out_scores

def process_goutput(*args, **kwargs
                   ) -> Tuple[Sample, TaxTree, SampleDataById,
                              SampleStats, Err]:
    """
    Process GO output files (to be usually called in parallel!).
    """
    # timing initialization
    start_time: float = time.perf_counter()
    # Recover input and parameters
    target_file: Filename = args[0]
    debug: bool = kwargs['debug']
    is_ctrl: bool = args[1]
    if debug:
        print(gray('Processing'), blue('ctrl' if is_ctrl else 'sample'),
              target_file, gray('...'))
        sys.stdout.flush()
    ontology: Ontology = kwargs['ontology']
    minhit: int = (kwargs['ctrlminhit'] if is_ctrl else kwargs['minhit'])
    minscore: Score = kwargs['ctrlminscore'] if is_ctrl else kwargs['minscore']
    including: Union[Tuple, Set[Id]] = ontology.including
    excluding: Union[Tuple, Set[Id]] = ontology.excluding
    scoring: Scoring = kwargs['scoring']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    sample: Sample = Sample(os.path.splitext(target_file)[0])
    error: Err = Err.NO_ERROR
    # Read GO output files to get abundances
    log: str
    counts: Counter[Id]
    scores: Dict[Id, Score]
    ranks: Ranks = Ranks({})
#    log, stat, counts, scores = read_blast_go(target_file, scoring, minscore)
    log, stat, counts, scores = read_kraken_output(
        target_file, scoring, minscore)
    output.write(log)
    # Update field in stat about control nature of the sample
    stat.is_ctrl = is_ctrl

    # Remove root counts, in case
    if kwargs['root'] and counts[ontology.ROOT]:
        vwrite(gray('Removing'), counts[ontology.ROOT],
               gray('"ROOT" reads... '))
        stat.seq = stat.seq._replace(filt=stat.seq.filt-counts[ontology.ROOT])
        stat.decrease_filtered_taxids()
        counts[ontology.ROOT] = 0
        scores[ontology.ROOT] = NO_SCORE
        vwrite(green('OK!'), '\n')

    # Building ontology tree
    output.write(gray('Building from raw data... '))
    vwrite(gray('\n  Building ontology tree with all-in-1... '))
    tree = TaxTree()
    ancestors: Set[Id]
    orphans: Set[Id]
    ancestors, orphans = ontology.get_ancestors(counts.keys())
    out = SampleDataById(['all'])
    tree.allin1(ontology=ontology, counts=counts, scores=scores,
                ancestors=ancestors, min_taxa=minhit,
                include=including, exclude=excluding, out=out)
    out.purge_counters()
    vwrite(green('OK!'), '\n')

    # Stats: Complete final value for TaxIDs after tree building and folding
    final_goids: int = len(out.counts) if out.counts is not None else 0
    stat.set_final_taxids(final_goids)

    # Check for additional loss of reads (due to include/exclude an orphans)
    output.write(gray('  Check for more seqs lost ([in/ex]clude affects)... '))
    if out.counts is not None:
        discard: int = sum(counts.values()) - sum(out.counts.values())
        if discard:
            output.write(blue('\n  Info:') + f' {discard} ' +
                         gray('additional seqs discarded (') +
                         f'{discard/sum(counts.values()):.3%} ' +
                         gray('of accepted)\n'))
        else:
            output.write(green('OK!\n'))
    else:
        output.write(red('No counts in sample tree!\n'))
    # Warn or give detailed stats about orphan GOid and orphan seqs
    if debug:
        vwrite(gray('  Checking GOid loss (orphans)... '))
        lost: int = 0
        if orphans:
            for orphan in orphans:
                vwrite(yellow('  Warning!'), gray('Orphan GOid'),
                       f'{orphan}\n')
                lost += counts[orphan]
            vwrite(yellow('  WARNING!'),
                   f'{len(orphans)} orphan GOids ('
                   f'{len(orphans)/len(counts):.2%} of accepted)\n'
                   f'    and {lost} orphan sequences ('
                   f'{lost/sum(counts.values()):.3%} of accepted)\n')
        else:
            vwrite(green('OK!\n'))
    elif orphans:
        output.write(yellow('\n  Warning!') + f' {len(orphans)} orphan GOids'
                     + gray(' (rerun with --debug for details)\n'))
    # Check the removal of TaxIDs (accumulation of leaves in parents)
    if debug and not excluding and including == {ontology.ROOT}:
        vwrite(gray('  Assess accumulation due to "folding the tree"...\n'))
        migrated: int = 0
        if out.counts is not None:
            for goid in counts:
                if out.counts[goid] == 0:
                    migrated += 1
                    vwrite(blue('  Info:'), gray(f'Folded GOid {goid} (') +
                           f'{ontology.get_name(goid)}' + gray(') with ') +
                           f'{counts[goid]}' + gray(' original seqs\n'))
        if migrated:
            vwrite(
                blue('  INFO:'), f'{migrated} GOids folded ('
                f'{migrated/len(+counts):.2%} of GOAF —GOids After Filtering—)'
                '\n')
            vwrite(
                blue('  INFO:'), f'Final assigned GOids: {final_goids} '
                f'(reduced to {final_goids/len(+counts):.2%} of '
                'number of GOAF)\n')
        else:
            vwrite(blue('  INFO:'), gray('No migration!'), green('OK!\n'))
    # Print last message and check if the sample is void
    if out.counts:
        output.write(sample + blue(' ctrl ' if is_ctrl else ' sample ')
                     + green('OK!\n'))
    elif is_ctrl:
        output.write(sample + red(' ctrl VOID!\n'))
        error = Err.VOID_CTRL
    else:
        output.write(sample + blue(' sample ') + yellow('VOID\n'))
        error = Err.VOID_SAMPLE

    # Timing results
    output.write(gray('Load elapsed time: ') +
                 f'{time.perf_counter() - start_time:.3g}' + gray(' sec\n'))
    print(output.getvalue())
    sys.stdout.flush()
    return sample, tree, out, stat, error

