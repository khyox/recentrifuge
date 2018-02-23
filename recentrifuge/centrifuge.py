"""
Functions directly related with Centrifuge (or Kraken).

"""

import collections as col
import io
import os
import sys
import time
from math import log10
from statistics import mean
from typing import Tuple, Counter, Callable, Optional, Set, Dict, List

from Bio import SeqIO

from recentrifuge.config import Filename, TaxId, Score, Scoring, Sample
from recentrifuge.config import UNCLASSIFIED, ROOT, NO_SCORE, Err, SampleStats
from recentrifuge.config import gray, red, green, yellow, blue
from recentrifuge.lmat import read_lmat_output
from recentrifuge.rank import Rank, Ranks
from recentrifuge.taxonomy import Taxonomy
from recentrifuge.trees import TaxTree, SampleDataByTaxId


def read_report(report_file: str) -> Tuple[str, Counter[TaxId],
                                           Dict[TaxId, Rank]]:
    """
    Read Centrifuge/Kraken report file

    Args:
        report_file: report file name

    Returns:
        log string, abundances counter, taxlevel dict

    """
    output: io.StringIO = io.StringIO(newline='')
    abundances: Counter[TaxId] = col.Counter()
    level_dic = {}
    output.write(f'\033[90mLoading report file {report_file}...\033[0m')
    try:
        with open(report_file, 'r') as file:
            for report_line in file:
                _, _, taxnum, taxlev, _tid, _ = report_line.split('\t')
                tid = TaxId(_tid)
                abundances[tid] = int(taxnum)
                level_dic[tid] = Rank.centrifuge(taxlev)
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        report_file + '"')
    else:
        output.write('\033[92m OK! \033[0m\n')
    return output.getvalue(), abundances, level_dic


def process_report(*args, **kwargs
                   ) -> Tuple[Sample, TaxTree, SampleDataByTaxId,
                              SampleStats, Err]:
    """
    Process Centrifuge/Kraken report files (to be usually called in parallel!).
    """
    # TODO: Full review to report support
    # Recover input and parameters
    filerep: Filename = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    mintaxa: int = kwargs['mintaxa']
    collapse: bool = taxonomy.collapse
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    debug: bool = kwargs['debug']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    sample: Sample = Sample(filerep)

    # Read Centrifuge/Kraken report file to get abundances
    log: str
    abundances: Counter[TaxId]
    log, abundances, _ = read_report(filerep)
    output.write(log)
    # Remove root counts, in case
    if kwargs['root']:
        vwrite(gray('Removing'), abundances[ROOT], gray('"ROOT" reads... '))
        abundances[ROOT] = 0
        vwrite(green('OK!'), '\n')

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy=taxonomy,
              counts=abundances)  # Grow tax tree from root node
    output.write('\033[92m OK! \033[0m\n')

    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    tree.prune(mintaxa, None, collapse, debug)
    tree.shape()
    output.write('\033[92m OK! \033[0m\n')

    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    new_abund: Counter[TaxId] = col.Counter()
    new_accs: Counter[TaxId] = col.Counter()
    ranks: Ranks = Ranks({})
    tree.get_taxa(abundance=new_abund,
                  accs=new_accs,
                  ranks=ranks,
                  mindepth=0, maxdepth=0,
                  include=including,
                  exclude=excluding)
    new_abund = +new_abund  # remove zero and negative counts
    if including or excluding:  # Recalculate accumulated counts
        new_tree = TaxTree()
        new_tree.grow(taxonomy, new_abund)  # Grow tree with new abund
        new_tree.shape()
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_tree.get_taxa(new_abund, new_accs)  # Get new accumulated counts
    out: SampleDataByTaxId = SampleDataByTaxId()
    out.set(counts=new_abund, ranks=ranks, accs=new_accs)
    output.write('\033[92m OK! \033[0m\n')
    print(output.getvalue())
    sys.stdout.flush()
    return sample, tree, out, SampleStats(), Err.NO_ERROR


def read_output(output_file: Filename,
                scoring: Scoring = Scoring.SHEL,
                minscore: Score = None,
                ) -> Tuple[str, SampleStats,
                           Counter[TaxId], Dict[TaxId, Score]]:
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
    all_scores: Dict[TaxId, List[Score]] = {}
    all_length: Dict[TaxId, List[int]] = {}
    num_read: int = 0
    nt_read: int = 0
    num_uncl: int = 0
    output.write(gray(f'Loading output file {output_file}... '))
    try:
        with open(output_file, 'r') as file:
            file.readline()  # discard header
            for output_line in file:
                _, _, _tid, _score, _, _, _length, *_ = output_line.split('\t')
                tid = TaxId(_tid)
                num_read += 1
                try:
                    # From Centrifuge score get "single hit equivalent length"
                    shel = Score(float(_score) ** 0.5 + 15)
                    length = int(_length)
                except ValueError:
                    print(red('Error'), f'parsing score ({_score}) for query',
                          f'length ({_length}) for taxid {_tid}',
                          f'in {output_file}. Ignoring line!')
                    continue
                nt_read += length
                if tid == UNCLASSIFIED:  # Just count unclassified reads
                    num_uncl += 1
                    continue
                elif minscore is not None and shel < minscore:
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
        raise Exception(red('\nERROR!'), f'Cannot read "{output_file}"')
    counts: Counter[TaxId] = Counter({tid: len(all_scores[tid])
                                      for tid in all_scores})
    output.write(green('OK!\n'))
    # Get statistics
    stat: SampleStats = SampleStats(
        minscore=minscore, nt_read=nt_read, scores=all_scores, lens=all_length,
        seq_read=num_read, seq_unclas=num_uncl,
        seq_filt=sum([len(scores) for scores in all_scores.values()]),
    )
    # Output statistics
    output.write(gray('  Seqs read: ') + f'{stat.seq.read:_d}\t' + gray('[')
                 + f'{stat.nt_read}' + gray(']\n'))
    output.write(gray('  Seqs clas: ') + f'{stat.seq.clas:_d}\t' + gray('(') +
                 f'{stat.get_unclas_ratio():.2%}' + gray(' unclassified)\n'))
    output.write(gray('  Seqs pass: ') + f'{stat.seq.filt:_d}\t' + gray('(') +
                 f'{stat.get_reject_ratio():.2%}' + gray(' rejected)\n'))
    output.write(gray('  Scores: min = ') + f'{stat.sco.mini:.1f},' +
                 gray(' max = ') + f'{stat.sco.maxi:.1f},' +
                 gray(' avr = ') + f'{stat.sco.mean:.1f}\n')

    output.write(gray('  Length: min = ') + f'{stat.len.mini},' +
                 gray(' max = ') + f'{stat.len.maxi},' +
                 gray(' avr = ') + f'{stat.len.mean}\n')
    output.write(f'  {stat.num_taxa}' + gray(f' taxa with assigned reads\n'))
    # Select score output
    out_scores: Dict[TaxId, Score]
    if scoring is Scoring.SHEL:
        out_scores = {tid: Score(mean(all_scores[tid])) for tid in all_scores}
    elif scoring is Scoring.LENGTH:
        out_scores = {tid: Score(mean(all_length[tid])) for tid in all_length}
    elif scoring is Scoring.LOGLENGTH:
        out_scores = {tid: Score(log10(mean(all_length[tid])))
                      for tid in all_length}
    elif scoring is Scoring.NORMA:
        scores: Dict[TaxId, Score] = {tid: Score(mean(all_scores[tid]))
                                      for tid in all_scores}
        lengths: Dict[TaxId, Score] = {tid: Score(mean(all_length[tid]))
                                       for tid in all_length}
        out_scores = {tid: Score(scores[tid] / lengths[tid] * 100)
                      for tid in scores}
    else:
        raise Exception(f'\n\033[91mERROR!\033[0m Unknown Scoring "{scoring}"')
    # Return
    return output.getvalue(), stat, counts, out_scores


def process_output(*args, **kwargs
                   ) -> Tuple[Sample, TaxTree, SampleDataByTaxId,
                              SampleStats, Err]:
    """
    Process Centrifuge/LMAT output files (to be usually called in parallel!).
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
    taxonomy: Taxonomy = kwargs['taxonomy']
    mintaxa: int = kwargs['ctrlmintaxa'] if is_ctrl else kwargs['mintaxa']
    minscore: Score = kwargs['ctrlminscore'] if is_ctrl else kwargs['minscore']
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    scoring: Scoring = kwargs['scoring']
    lmat: bool = kwargs['lmat']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    sample: Sample = Sample(os.path.splitext(target_file)[0])
    error: Err = Err.NO_ERROR
    # Read Centrifuge/LMAT output files to get abundances
    read_method: Callable[
        [Filename, Scoring, Optional[Score]],  # Input
        Tuple[str, SampleStats, Counter[TaxId], Dict[TaxId, Score]]  # Output
    ]
    if lmat:
        read_method = read_lmat_output
    else:
        read_method = read_output
    log: str
    counts: Counter[TaxId]
    scores: Dict[TaxId, Score]
    log, stat, counts, scores = read_method(target_file, scoring, minscore)
    output.write(log)
    # Remove root counts, in case
    if kwargs['root']:
        vwrite(gray('Removing'), counts[ROOT], gray('"ROOT" reads... '))
        counts[ROOT] = 0
        scores[ROOT] = NO_SCORE
        vwrite(green('OK!'), '\n')

    # Building taxonomy tree
    output.write(gray('Building from raw data... '))
    vwrite(gray('\n  Building taxonomy tree with all-in-1... '))
    tree = TaxTree()
    ancestors: Set[TaxId]
    orphans: Set[TaxId]
    ancestors, orphans = taxonomy.get_ancestors(counts.keys())
    out = SampleDataByTaxId(['all'])
    tree.allin1(taxonomy=taxonomy, counts=counts, scores=scores,
                ancestors=ancestors, min_taxa=mintaxa,
                include=including, exclude=excluding, out=out)
    out.purge_counters()
    vwrite(green('OK!'), '\n')

    # Give stats about orphan taxid
    if debug:
        vwrite(gray('  Checking taxid loss (orphans)... '))
        lost: int = 0
        if orphans:
            for orphan in orphans:
                vwrite(yellow('Warning!'), f'Orphan taxid={orphan}\n')
                lost += counts[orphan]
            vwrite(yellow('WARNING!'),
                   f'{len(orphans)} orphan taxids ('
                   f'{len(orphans)/len(counts):.2%} of total)\n'
                   f'{lost} orphan sequences ('
                   f'{lost/sum(counts.values()):.3%} of total)\n')
        else:
            vwrite(green('OK!\n'))
    # Check the lost of taxids (plasmids typically) under some conditions
    if debug and not excluding and not including:
        vwrite(gray('  Additional checking of taxid loss... '))
        lost = 0
        for taxid in counts:
            if not out.counts[taxid]:
                lost += 1
                vwrite(yellow('Warning!'), f'Lost taxid={taxid}: '
                                           f'{taxonomy.get_name(taxid)}\n')
        if lost:
            vwrite(yellow('WARNING!'), f'Lost {lost} taxids ('
                                       f'{lost/len(counts):.2%} of total)'
                                       '\n')
        else:
            vwrite(green('OK!\n'))

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
