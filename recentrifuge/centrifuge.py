"""
Functions directly related with Centrifuge (or Kraken).

"""

import collections as col
import io
import sys
from typing import Tuple, Counter, Set, Dict

from recentrifuge.config import ROOT, Filename, TaxId
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
                              Counter[TaxId], Counter[TaxId]]:
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
    log, abundances, _ = read_report(filerep)
    output.write(log)

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy, abundances, ROOT, [])  # Grow tax tree from root node
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
    tree.get_taxa(new_abund, new_accs, ranks,
                  mindepth=0, maxdepth=0,
                  include=including,  # ('2759',),  # ('135613',),
                  exclude=excluding)  # ('9606',))  # ('255526',))
    new_abund = +new_abund  # remove zero and negative counts
    if including or excluding:  # Recalculate accumulated counts
        new_tree = TaxTree()
        new_tree.grow(taxonomy, new_abund, ROOT, [])  # Grow tree with new abund
        new_tree.accumulate()
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_tree.get_taxa(new_abund, new_accs)  # Get new accumulated counts
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    output.write('\033[92m OK! \033[0m\n')
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return sample, tree, taxlevels, new_abund, new_accs