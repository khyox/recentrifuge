"""
Common functions to taxonomic classifiers.

"""

import io
import os
import sys
import time
import collections as col
from typing import Tuple, Set, Callable, Optional, Counter, Dict

from recentrifuge.centrifuge import read_output, read_report
from recentrifuge.clark import read_clark_output
from recentrifuge.config import Sample, Err, Filename, Score, Id
from recentrifuge.stats import SampleStats
from recentrifuge.config import Scoring, Classifier
from recentrifuge.config import CELLULAR_ORGANISMS, NO_SCORE
from recentrifuge.config import gray, blue, green, yellow, red
from recentrifuge.lmat import read_lmat_output
from recentrifuge.ontology import Ontology
from recentrifuge.rank import Ranks
from recentrifuge.trees import TaxTree, SampleDataById


def process_output(*args, **kwargs
                   ) -> Tuple[Sample, TaxTree, SampleDataById,
                              SampleStats, Err]:
    """
    Process classifiers output files (to be usually called in parallel!).
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
    mintaxa: int = kwargs['ctrlmintaxa'] if is_ctrl else kwargs['mintaxa']
    minscore: Score = kwargs['ctrlminscore'] if is_ctrl else kwargs['minscore']
    including: Set[Id] = ontology.including
    excluding: Set[Id] = ontology.excluding
    scoring: Scoring = kwargs['scoring']
    classifier: Classifier = kwargs['classifier']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    sample: Sample = Sample(os.path.splitext(target_file)[0])
    error: Err = Err.NO_ERROR
    # Read taxonomic classifier output files to get abundances
    read_method: Callable[
        [Filename, Scoring, Optional[Score]],  # Input
        Tuple[str, SampleStats, Counter[Id], Dict[Id, Score]]  # Output
    ]
    if classifier is Classifier.CLARK:
        read_method = read_clark_output
    elif classifier is Classifier.LMAT:
        read_method = read_lmat_output
    elif classifier is Classifier.CENTRIFUGE:
        read_method = read_output
    else:
        raise Exception(f'taxclass: Unknown taxon classifier "{classifier}".')
    log: str
    counts: Counter[Id]
    scores: Dict[Id, Score]
    log, stat, counts, scores = read_method(target_file, scoring, minscore)
    output.write(log)
    # Update field in stat about control nature of the sample
    stat.is_ctrl = is_ctrl
    # Move cellular_organisms counts to root, in case
    if ontology.collapse and counts[CELLULAR_ORGANISMS]:
        vwrite(gray('Moving'), counts[CELLULAR_ORGANISMS],
               gray('"CELLULAR_ORGANISMS" reads to "ROOT"... '))
        if counts[ontology.ROOT]:
            stat.num_taxa -= 1
            scores[ontology.ROOT] = (
                    (scores[CELLULAR_ORGANISMS] * counts[CELLULAR_ORGANISMS] +
                     scores[ontology.ROOT] * counts[ontology.ROOT])
                    / (counts[CELLULAR_ORGANISMS] + counts[ontology.ROOT]))
        else:
            scores[ontology.ROOT] = scores[CELLULAR_ORGANISMS]
        counts[ontology.ROOT] += counts[CELLULAR_ORGANISMS]
        counts[CELLULAR_ORGANISMS] = 0
        scores[CELLULAR_ORGANISMS] = NO_SCORE
    # Remove root counts, in case
    if kwargs['root'] and counts[ontology.ROOT]:
        vwrite(gray('Removing'), counts[ontology.ROOT],
               gray('"ROOT" reads... '))
        stat.seq = stat.seq._replace(filt=stat.seq.filt-counts[ontology.ROOT])
        stat.num_taxa -= 1
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
                                           f'{ontology.get_name(taxid)}\n')
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


def process_report(*args, **kwargs
                   ) -> Tuple[Sample, TaxTree, SampleDataById,
                              SampleStats, Err]:
    """
    Process Centrifuge/Kraken report files (to be usually called in parallel!).
    """
    # TODO: Full review to report support
    # Recover input and parameters
    filerep: Filename = args[0]
    ontology: Ontology = kwargs['ontology']
    mintaxa: int = kwargs['mintaxa']
    collapse: bool = ontology.collapse
    including: Set[Id] = ontology.including
    excluding: Set[Id] = ontology.excluding
    debug: bool = kwargs['debug']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args):
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    sample: Sample = Sample(filerep)

    # Read Centrifuge/Kraken report file to get abundances
    log: str
    abundances: Counter[Id]
    log, abundances, _ = read_report(filerep)
    output.write(log)
    # Remove root counts, in case
    if kwargs['root']:
        vwrite(gray('Removing'), abundances[ontology.ROOT],
               gray('"ROOT" reads... '))
        abundances[ontology.ROOT] = 0
        vwrite(green('OK!'), '\n')

    # Build ontology tree
    output.write('  \033[90mBuilding ontology tree...\033[0m')
    tree = TaxTree()
    tree.grow(ontology=ontology,
              counts=abundances)  # Grow tax tree from root node
    output.write('\033[92m OK! \033[0m\n')

    # Prune the tree
    output.write('  \033[90mPruning ontology tree...\033[0m')
    tree.prune(mintaxa, None, collapse, debug)
    tree.shape()
    output.write('\033[92m OK! \033[0m\n')

    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    new_abund: Counter[Id] = col.Counter()
    new_accs: Counter[Id] = col.Counter()
    ranks: Ranks = Ranks({})
    tree.get_taxa(counts=new_abund,
                  accs=new_accs,
                  ranks=ranks,
                  mindepth=0, maxdepth=0,
                  include=including,
                  exclude=excluding)
    new_abund = +new_abund  # remove zero and negative counts
    if including or excluding:  # Recalculate accumulated counts
        new_tree = TaxTree()
        new_tree.grow(ontology, new_abund)  # Grow tree with new abund
        new_tree.shape()
        new_abund = col.Counter()  # Reset abundances
        new_accs = col.Counter()  # Reset accumulated
        new_tree.get_taxa(new_abund, new_accs)  # Get new accumulated counts
    out: SampleDataById = SampleDataById()
    out.set(counts=new_abund, ranks=ranks, accs=new_accs)
    output.write('\033[92m OK! \033[0m\n')
    print(output.getvalue())
    sys.stdout.flush()
    return sample, tree, out, SampleStats(), Err.NO_ERROR
