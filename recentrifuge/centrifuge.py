"""
Functions directly related with Centrifuge (or Kraken).

"""

import collections as col
import io
import sys
from math import log10
from statistics import mean
from typing import Tuple, Counter, Set, Dict, List

from Bio import SeqIO

from recentrifuge.LMAT import read_lmat_output
from recentrifuge.config import Filename, TaxId, Score, UNCLASSIFIED, Scoring
from recentrifuge.config import nucleotides, gray, red, green
from recentrifuge.core import Rank, Taxonomy, TaxTree, Ranks, TaxLevels, Sample


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


def process_report(*args,
                   **kwargs
                   ) -> Tuple[Sample, TaxTree, TaxLevels,
                              Counter[TaxId], Counter[TaxId], None]:
    """
    Process Centrifuge/Kraken report files (to be usually called in parallel!).
    """
    # Recover input and parameters
    filerep: Filename = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    mintaxa: int = kwargs['mintaxa']
    collapse: bool = taxonomy.collapse
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    debug: bool = kwargs['debug']
    output: io.StringIO = io.StringIO(newline='')

    sample: Sample = Sample(filerep)

    # Read Centrifuge/Kraken report file to get abundances
    log: str
    abundances: Counter[TaxId]
    log, abundances, _ = read_report(filerep)
    output.write(log)

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy=taxonomy,
              abundances=abundances)  # Grow tax tree from root node
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
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    output.write('\033[92m OK! \033[0m\n')
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return sample, tree, taxlevels, new_abund, new_accs, None


def read_output(output_file: Filename,
                scoring: Scoring = Scoring.SHEL,
                minscore: Score = None,
                ) -> Tuple[str, Counter[TaxId],
                           Dict[TaxId, Score]]:
    """
    Read Centrifuge output file

    Args:
        output_file: output file name
        scoring: type of scoring to be applied (see Scoring class)
        minscore: minimum confidence level for the classification

    Returns:
        log string, abundances counter, scores dict

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
                elif minscore is not None:
                    if shel < minscore:  # Ignore read if low confidence
                        continue
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
    abundances: Counter[TaxId] = Counter({tid: len(all_scores[tid])
                                          for tid in all_scores})
    output.write(green('OK!\n'))
    # Basic output statistics
    num_clas: int = num_read - num_uncl
    num_seqs: int = sum([len(scores) for scores in all_scores.values()])
    output.write(gray('  Seqs read: ') + f'{num_read:_d}\t' + gray('[') +
                 nucleotides(nt_read) + gray(']\n'))
    output.write(gray('  Seqs clas: ') + f'{num_clas:_d}\t' + gray('(') +
                 f'{num_uncl/num_read:.2%}' + gray(' unclassified)\n'))
    output.write(gray('  Seqs pass: ') + f'{num_seqs:_d}\t' + gray('(') +
                 f'{1-(num_seqs/num_clas):.2%}' + gray(' filtered)\n'))
    max_score = max([max(s) for s in all_scores.values()])
    min_score = min([min(s) for s in all_scores.values()])
    mean_score = mean([mean(s) for s in all_scores.values()])
    output.write(f'\033[90m  Scores: min =\033[0m {min_score:.1f},'
                 f'\033[90m max =\033[0m {max_score:.1f},'
                 f'\033[90m avr =\033[0m {mean_score:.1f}\n')
    max_length = max([max(l) for l in all_length.values()])
    min_length = min([min(l) for l in all_length.values()])
    mean_length = mean([mean(l) for l in all_length.values()])
    output.write(f'\033[90m  Length: min =\033[0m {int(min_length)},'
                 f'\033[90m max =\033[0m {int(max_length)},'
                 f'\033[90m avr =\033[0m {mean_length:.1f}\n')
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
    return output.getvalue(), abundances, out_scores


def process_output(*args,
                   **kwargs
                   ) -> Tuple[Sample, TaxTree, TaxLevels,
                              Counter[TaxId], Counter[TaxId],
                              Dict[TaxId, Score]]:
    """
    Process Centrifuge/LMAT output files (to be usually called in parallel!).
    """
    # Recover input and parameters
    fileout: Filename = args[0]
    is_ctrl: bool = args[1]
    taxonomy: Taxonomy = kwargs['taxonomy']
    mintaxa: int = kwargs['ctrlmintaxa'] if is_ctrl else kwargs['mintaxa']
    minscore: Score = kwargs['ctrlminscore'] if is_ctrl else kwargs['minscore']
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    debug: bool = kwargs['debug']
    scoring: Scoring = kwargs['scoring']
    lmat: bool = kwargs['lmat']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*a, **k):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(''.join(*a, **k))

    sample: Sample = Sample(fileout)

    # Read Centrifuge/LMAT output files to get abundances
    if lmat:
        read_method = read_lmat_output
    else:
        read_method = read_output
    log: str
    abundances: Counter[TaxId]
    scores: Dict[TaxId, Score]
    log, abundances, scores = read_method(fileout, scoring, minscore)
    output.write(log)

    # Build taxonomy tree
    vwrite(gray('  Building taxonomy tree...'))
    tree = TaxTree()
    tree.grow(taxonomy=taxonomy,
              abundances=abundances,
              scores=scores)  # Grow tax tree from root node
    tree.shape()
    vwrite(green('OK!'), '\n')

    if debug:  # Check the lost of taxids (plasmids typically)
        output.write('  \033[90mChecking taxid loss...\033[0m')
        n_abund: Counter[TaxId] = col.Counter()
        n_accs: Counter[TaxId] = col.Counter()
        n_scores: Dict[TaxId, Score] = {}
        n_ranks: Ranks = Ranks({})
        tree.get_taxa(n_abund, n_accs, n_scores, n_ranks,
                      mindepth=0, maxdepth=0)
        lost: int = 0
        for taxid in abundances:
            if not n_abund[taxid]:
                lost += 1
                output.writelines(
                    f'\033[93mWarning!\033[0m Lost '
                    f'taxid={taxid}: {taxonomy.get_name(taxid)}\n')
        if lost:
            output.write(f'\033[93mWARNING!\033[0m Lost {lost} taxids'
                         f' ({lost/len(abundances):.2%} of total)\n')
        else:
            output.write(green('OK!\n'))

    # Prune the tree
    vwrite(gray('  Pruning taxonomy tree...'))
    tree.prune(min_taxa=mintaxa, debug=debug)
    vwrite(green('OK!'), '\n')

    # Get the taxa with their abundances, scores and taxonomical levels
    vwrite(gray('  Filtering taxa...'))
    new_abund: Counter[TaxId] = col.Counter()
    new_accs: Counter[TaxId] = col.Counter()
    new_scores: Dict[TaxId, Score] = {}
    ranks: Ranks = Ranks({})
    tree.get_taxa(new_abund, new_accs, new_scores, ranks,
                  mindepth=0, maxdepth=0,
                  include=including,  # ('2759',),  # ('135613',),
                  exclude=excluding)  # ('9606',))  # ('255526',))
    new_abund = +new_abund  # remove zero and negative counts
    if including or excluding:  # Recalculate accumulated counts
        new_tree = TaxTree()
        new_tree.grow(taxonomy=taxonomy,
                      abundances=new_abund,
                      scores=new_scores)  # Grow tree with new abund
        new_tree.shape()
        new_tree.prune(min_taxa=mintaxa, debug=debug)
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_scores = {}  # Reset score
        new_tree.get_taxa(new_abund, new_accs, new_scores)  # Get new values
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    vwrite(green('OK!'), '\n')
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return sample, tree, taxlevels, new_abund, new_accs, new_scores


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
