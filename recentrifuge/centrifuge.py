"""
Functions directly related with Centrifuge (or Kraken).

"""

import collections as col
import io
import sys
from statistics import mean
from typing import Tuple, Counter, Set, Dict, List

from Bio import SeqIO

from recentrifuge.config import Filename, TaxId, Score, UNCLASSIFIED
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
    verb: bool = kwargs['verb']
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
    tree.prune(mintaxa, None, collapse, verb)
    tree.accumulate()
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
        new_tree.accumulate()
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_tree.get_taxa(new_abund, new_accs)  # Get new accumulated counts
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    output.write('\033[92m OK! \033[0m\n')
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return sample, tree, taxlevels, new_abund, new_accs, None


def read_output(output_file: str) -> Tuple[str, Counter[TaxId],
                                           Dict[TaxId, Score]]:
    """
    Read Centrifuge output file

    Args:
        output_file: output file name

    Returns:
        log string, abundances counter, scores dict

    """
    output: io.StringIO = io.StringIO(newline='')
    all_scores: Dict[TaxId, List[Score]] = {}
    output.write(f'\033[90mLoading output file {output_file}...\033[0m')
    try:
        with open(output_file, 'r') as file:
            file.readline()  # discard header
            for output_line in file:
                _, _, _tid, _score, *_ = output_line.split('\t')
                tid = TaxId(_tid)
                # From Centrifuge score get the "single hit equivalent length"
                try:
                    score = Score(float(_score)**0.5 + 15)
                except ValueError:
                    print(f'Error parsing score ({_score}) for taxid {_tid}'
                          f' in file {output_file}...')
                    raise
                try:
                    all_scores[tid].append(score)
                except KeyError:
                    all_scores[tid] = [score, ]
    except FileNotFoundError:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        output_file + '"')
    abundances: Counter[TaxId] = Counter({tid: len(all_scores[tid])
                                          for tid in all_scores})
    scores: Dict[TaxId, Score] = {tid: Score(mean(all_scores[tid]))
                                  for tid in all_scores}
    output.write('\033[92m OK! \033[0m\n')
    # Basic output statistics
    num_seqs: int = sum([len(scores) for scores in all_scores.values()])
    num_uncl: int = len(all_scores.pop(UNCLASSIFIED, []))
    output.write(f'  \033[90mSeqs read: \033[0m{num_seqs:_d} \033[90m\t'
                 f'(\033[0m{num_uncl/num_seqs:.1%}\033[90m unclassified)\n')
    max_score = max([max(scores) for scores in all_scores.values()])
    min_score = min([min(scores) for scores in all_scores.values()])
    mean_score = mean([mean(scores) for scores in all_scores.values()])
    output.write(f'\033[90m  Scores: min =\033[0m {int(min_score)},'
                 f'\033[90m max =\033[0m {int(max_score)},'
                 f'\033[90m avr =\033[0m {int(mean_score)}\n')

    return output.getvalue(), abundances, scores


def process_output(*args,
                   **kwargs
                   ) -> Tuple[Sample, TaxTree, TaxLevels,
                              Counter[TaxId], Counter[TaxId],
                              Dict[TaxId, Score]]:
    """
    Process Centrifuge output files (to be usually called in parallel!).
    """
    # Recover input and parameters
    fileout: Filename = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    mintaxa: int = kwargs['mintaxa']
    collapse: bool = taxonomy.collapse
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    verb: bool = kwargs['verb']
    output: io.StringIO = io.StringIO(newline='')

    sample: Sample = Sample(fileout)

    # Read Centrifuge/Kraken report file to get abundances
    log: str
    abundances: Counter[TaxId]
    scores: Dict[TaxId, Score]
    log, abundances, scores = read_output(fileout)
    output.write(log)

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy=taxonomy,
              abundances=abundances,
              scores=scores)  # Grow tax tree from root node
    output.write('\033[92m OK! \033[0m\n')

    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    tree.prune(mintaxa=mintaxa, verb=verb)
    tree.accumulate()
    output.write('\033[92m OK! \033[0m\n')

    # Get the taxa with their abundances, scores and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
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
        new_tree.accumulate()
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_scores = {}  # Reset score
        new_tree.get_taxa(new_abund, new_accs, new_scores)  # Get new values
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    output.write('\033[92m OK! \033[0m\n')
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
