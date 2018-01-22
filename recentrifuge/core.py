"""
Recentrifuge core functions.

"""
import collections as col
import copy
import csv
import io
import statistics
import subprocess
import sys
from typing import List, Set, Counter, Tuple, Union, Dict

from recentrifuge.config import Filename, Sample, TaxId
from recentrifuge.config import HTML_SUFFIX, CELLULAR_ORGANISMS, ROOT, EPS
from recentrifuge.config import Parents, Score
from recentrifuge.config import STR_CONTROL, STR_EXCLUSIVE, STR_SHARED
from recentrifuge.config import STR_SHARED_CONTROL
from recentrifuge.rank import Rank, TaxLevels
from recentrifuge.shared_counter import SharedCounter
from recentrifuge.taxonomy import Taxonomy
from recentrifuge.trees import TaxTree


def process_rank(*args,
                 **kwargs
                 ) -> Tuple[List[Sample],
                            Dict[Sample, Counter[TaxId]],
                            Dict[Sample, Counter[TaxId]],
                            Dict[Sample, Dict[TaxId, Score]]]:
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """

    # Recover input and parameters
    rank: Rank = args[0]
    controls: int = kwargs['controls']
    mintaxa = kwargs['mintaxa']
    taxonomy: Taxonomy = kwargs['taxonomy']
    including = taxonomy.including
    excluding = taxonomy.excluding
    trees: Dict[Sample, TaxTree] = kwargs['trees']
    taxids: Dict[Sample, TaxLevels] = kwargs['taxids']
    accs: Dict[Sample, Counter[TaxId]] = kwargs['accs']
    files: List[Sample] = kwargs['files']

    # Declare/define variables
    samples: List[Sample] = []
    abundances: Dict[Sample, Counter[TaxId]]
    scores: Dict[Sample, Dict[TaxId, Score]] = {}
    # pylint: disable = unused-variable
    shared_abundance: SharedCounter = SharedCounter()
    shared_score: SharedCounter = SharedCounter()
    shared_ctrl_abundance: SharedCounter = SharedCounter()
    shared_ctrl_score: SharedCounter = SharedCounter()
    # pylint: enable = unused-variable
    output: io.StringIO = io.StringIO(newline='')

    output.write(f'\033[90mAnalysis for taxonomic rank "'
                 f'\033[95m{rank.name.lower()}\033[90m":\033[0m\n')

    def cross_analysis(iteration, file):
        """Cross analysis: exclusive and part of shared&ctrl"""
        nonlocal shared_abundance, shared_score
        nonlocal shared_ctrl_abundance, shared_ctrl_score

        def generate_exclusive():
            """Generate a new tree to recalculate accs (and scores)"""
            nonlocal exclusive_score, sample
            tree = TaxTree()
            tree.grow(taxonomy=taxonomy,
                      abundances=exclusive_abundance,
                      scores=exclusive_score)
            tree.shape()
            exclusive_acc: Counter[TaxId] = col.Counter()
            exclusive_score = {}
            tree.get_taxa(accs=exclusive_acc,
                          scores=exclusive_score,
                          include=including,
                          exclude=excluding)
            sample = Sample(f'{STR_EXCLUSIVE}_{rank.name.lower()}_{file}')
            samples.append(sample)
            abundances[sample] = exclusive_abundance
            accs[sample] = exclusive_acc
            scores[sample] = exclusive_score

        def partial_shared_update(i):
            """Perform shared and shared-control taxa partial evaluations"""
            nonlocal shared_abundance, shared_score
            nonlocal shared_ctrl_abundance, shared_ctrl_score
            if i == 0:  # 1st iteration: Initialize shared abundance and score
                shared_abundance.update(sub_shared_abundance)
                shared_score.update(sub_shared_score)
            elif i < controls:  # Just update shared abundance and score
                shared_abundance &= sub_shared_abundance
                shared_score &= sub_shared_score
            elif i == controls:  # Initialize shared-control counters
                shared_abundance &= sub_shared_abundance
                shared_score &= sub_shared_score
                shared_ctrl_abundance.update(sub_shared_abundance)
                shared_ctrl_score.update(sub_shared_score)
            elif controls:  # Both: Accumulate shared abundance and score
                shared_abundance &= sub_shared_abundance
                shared_score &= sub_shared_score
                shared_ctrl_abundance &= sub_shared_abundance
                shared_ctrl_score &= sub_shared_score
            else:  # Both: Accumulate shared abundance and score (no controls)
                shared_abundance &= sub_shared_abundance
                shared_score &= sub_shared_score

        exclude: Set[TaxId] = set()
        # Get taxids at this rank that are present in the other samples
        for sample in (_file for _file in files if _file != file):
            exclude.update(taxids[sample][rank])
        exclude.update(excluding)  # Add explicit excluding taxa if any
        output.write(f'  \033[90mExclusive: From \033[0m{file}\033[90m '
                     f'excluding {len(exclude)} taxa. '
                     f'Generating sample...\033[0m')
        leveled_tree = copy.deepcopy(trees[file])
        leveled_tree.prune(mintaxa, rank)  # Prune the copy to the desired rank
        leveled_tree.shape()
        # Get abundance for the exclusive analysis
        exclusive_abundance: Counter[TaxId] = col.Counter()
        exclusive_score: Dict[TaxId, Score] = {}
        leveled_tree.get_taxa(abundance=exclusive_abundance,
                              scores=exclusive_score,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=exclude,  # Extended exclusion set
                              just_level=rank,
                              )
        exclusive_abundance = +exclusive_abundance  # remove counts <= 0
        if exclusive_abundance:  # Avoid adding empty samples
            generate_exclusive()
            output.write('\033[92m OK! \033[0m\n')
        else:
            output.write('\033[93m VOID \033[0m\n')
        # Get partial abundance and score for the shared analysis
        sub_shared_abundance: SharedCounter = SharedCounter()
        sub_shared_score: SharedCounter = SharedCounter()
        leveled_tree.get_taxa(abundance=sub_shared_abundance,
                              scores=sub_shared_score,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=excluding,
                              just_level=rank,
                              )
        # Scale scores by abundance
        sub_shared_score *= sub_shared_abundance
        partial_shared_update(iteration)

    def shared_analysis():
        """Perform last steps of shared taxa analysis"""
        nonlocal shared_abundance, shared_score
        # Normalize scaled scores by total abundance (after eliminating zeros)
        shared_score /= (+shared_abundance)
        # Get averaged abundance by number of samples
        shared_abundance //= len(files)
        tree = TaxTree()
        tree.grow(taxonomy=taxonomy,
                  abundances=shared_abundance,
                  scores=shared_score)
        tree.shape()
        new_abund: Counter[TaxId] = col.Counter()
        new_accs: Counter[TaxId] = col.Counter()
        new_score: Dict[TaxId, Score] = {}
        tree.get_taxa(new_abund, new_accs, new_score,
                      mindepth=0, maxdepth=0,
                      include=including,
                      exclude=excluding)
        new_abund = +new_abund
        output.write(f'  \033[90mShared: Including {len(new_abund)}'
                     f' shared taxa. Generating sample...\033[0m')
        sample = Sample(f'{STR_SHARED}_{rank.name.lower()}')
        samples.append(sample)
        abundances[Sample(sample)] = new_abund
        accs[Sample(sample)] = new_accs
        scores[sample] = new_score
        output.write('\033[92m OK! \033[0m\n')

    def control_analysis():
        """Perform last steps of control and shared controls analysis"""
        nonlocal shared_ctrl_abundance, shared_ctrl_score

        def generate_control(file_name):
            """Generate a new tree to be able to calculate accs and scores"""
            nonlocal tree, sample, control_score
            tree = TaxTree()
            tree.grow(taxonomy=taxonomy,
                      abundances=control_abundance,
                      scores=control_score)
            tree.shape()
            control_acc: Counter[TaxId] = col.Counter()
            control_score = {}
            tree.get_taxa(accs=control_acc,
                          scores=control_score,
                          include=including,
                          exclude=excluding)
            sample = Sample(f'{STR_CONTROL}_{rank.name.lower()}_{file_name}')
            samples.append(sample)
            abundances[sample] = control_abundance
            accs[sample] = control_acc
            scores[sample] = control_score

        def advanced_control_removal():
            """Implement advanced control removal"""
            nonlocal exclude_sets
            abund: List[float]
            abund_norm: List[float]
            crossover: List[bool]
            exclude_sets = {file: set() for file in files[controls:]}
            for tid in exclude_candidates:
                # Severe contaminant check
                is_severe_contaminant: bool = False
                for ctrl in files[:controls]:
                    # Severe contaminant if its rel.freq >1% in any control
                    assert accs[ctrl][ROOT] > 0, 'Panic! Empty tree!'
                    if accs[ctrl][tid] / accs[ctrl][ROOT] > 0.01:
                        is_severe_contaminant = True
                        for exclude_set in exclude_sets.values():
                            exclude_set.add(tid)
                        break
                if is_severe_contaminant:
                    continue
                # Broad/crossover contaminant checks
                abund = [accs[file][tid] for file in files[controls:]]
                mdn: float = statistics.median(abund)
                mad: float = statistics.median([abs(mdn - a) for a in abund])
                if mdn < EPS:
                    print('only-ctrl', tid, taxonomy.get_name(tid),
                          [accs[ctrl][tid] for ctrl in files[:controls]],
                          abund)
                elif mad < EPS:
                    print('no mad', tid, taxonomy.get_name(tid),
                          [accs[ctrl][tid] for ctrl in files[:controls]],
                          abund)
                else:
                    abund_norm = [(a - mdn) / mad for a in abund]
                    crossover = [(a > 2*mad) for a in abund_norm]
                    if any(crossover):
                        print('crossover', tid, taxonomy.get_name(tid),
                              abund_norm, crossover)
                    else:
                        print('contamination', tid, taxonomy.get_name(tid),
                              abund_norm, crossover)

        # Get taxids at this rank that are present in the control samples
        exclude_candidates = set()
        for i in range(controls):
            exclude_candidates.update(taxids[files[i]][rank])
        exclude_sets: Dict[Sample, Set[TaxId]]
        if len(samples) - controls > 0:
            advanced_control_removal()
        else:  # If less than 3 samples, just apply strict control
            exclude_sets = {file: exclude_candidates
                            for file in files[controls::]}
        # Add explicit excluding taxa (if any) to exclude sets
        for exclude_set in exclude_sets.values():
            exclude_set.update(excluding)
        exclude_candidates.update(excluding)
        # Process each sample excluding control taxa
        for file in files[controls:]:
            output.write(f'  \033[90mCtrl: From \033[0m{file}\033[90m '
                         f'excluding {len(exclude_sets[file])} ctrl taxa. '
                         f'Generating sample...\033[0m')
            leveled_tree = copy.deepcopy(trees[file])
            leveled_tree.prune(mintaxa, rank)
            leveled_tree.shape()
            control_abundance: Counter[TaxId] = col.Counter()
            control_score: Dict[TaxId, Score] = {}
            leveled_tree.get_taxa(abundance=control_abundance,
                                  scores=control_score,
                                  include=including,
                                  exclude=exclude_sets[file],
                                  just_level=rank,
                                  )
            control_abundance = +control_abundance  # remove counts <= 0
            if control_abundance:  # Avoid adding empty samples
                generate_control(file)
                output.write('\033[92m OK! \033[0m\n')
            else:
                output.write('\033[93m VOID \033[0m\n')
        # Shared-control taxa final analysis
        output.write('  \033[90mShared-ctrl:\033[0m ')
        if shared_ctrl_abundance:
            # Normalize scaled scores by total abundance
            shared_ctrl_score /= (+shared_ctrl_abundance)
            # Get averaged abundance by number of samples minus ctrl samples
            shared_ctrl_abundance //= (len(files) - controls)
            tree = TaxTree()
            tree.grow(taxonomy=taxonomy,
                      abundances=shared_ctrl_abundance,
                      scores=shared_ctrl_score)
            tree.shape()
            new_shared_ctrl_abundance: Counter[TaxId] = col.Counter()
            new_shared_ctrl_score: Dict[TaxId, Score] = {}
            tree.get_taxa(abundance=new_shared_ctrl_abundance,
                          scores=new_shared_ctrl_score,
                          include=including,
                          exclude=exclude_candidates,  # Extended exclusion set
                          just_level=rank,
                          )
            new_shared_ctrl_abundance = +new_shared_ctrl_abundance
            if new_shared_ctrl_abundance:  # Avoid adding empty samples
                output.write(
                    f'\033[90mIncluding '
                    f'{len(new_shared_ctrl_abundance)} shared taxa. '
                    f'Generating sample...\033[0m')

                # Generate a new tree to be able to calculate accs and scores
                tree = TaxTree()
                tree.grow(taxonomy=taxonomy,
                          abundances=new_shared_ctrl_abundance,
                          scores=new_shared_ctrl_score)
                tree.shape()
                new_shared_ctrl_acc: Counter[TaxId] = col.Counter()
                new_shared_ctrl_score = {}
                tree.get_taxa(accs=new_shared_ctrl_acc,
                              scores=new_shared_ctrl_score,
                              include=including,
                              exclude=excluding)
                sample = Sample(f'{STR_SHARED_CONTROL}_{rank.name.lower()}')
                samples.append(sample)
                abundances[sample] = new_shared_ctrl_abundance
                accs[sample] = new_shared_ctrl_acc
                scores[sample] = new_shared_ctrl_score
                output.write('\033[92m OK! \033[0m\n')
            else:
                output.write(f'\033[90mNo final shared-ctrl taxa!'
                             f'\033[93m VOID\033[90m sample.\033[0m \n')
        else:
            output.write(f'\033[90mNo shared-ctrl taxa!'
                         f'\033[93m VOID\033[90m sample.\033[0m \n')

    # Cross analysis iterating by output: exclusive and part of shared&ctrl
    for numfile, filename in zip(range(len(files)), files):
        cross_analysis(numfile, filename)

    # Shared taxa final analysis
    shared_abundance = +shared_abundance  # remove counts <= 0
    if shared_abundance:
        shared_analysis()
    else:
        output.write(f'  \033[90mShared: No shared taxa!'
                     f'\033[93m VOID\033[90m sample.\033[0m \n')

    # Control sample subtraction
    if controls:
        control_analysis()

    # Print output and return
    print(output.getvalue())
    sys.stdout.flush()
    return samples, abundances, accs, scores


def write_lineage(parents: Parents,
                  names: Dict[TaxId, str],
                  tree: TaxTree,
                  lineage_file: str,
                  nodes: Counter[TaxId],
                  collapse: bool = True,
                  ) -> str:
    """
    Writes a lineage file understandable by Krona.

    Args:
        parents: dictionary of parents for every TaxId.
        names: dictionary of names for every TaxId.
        tree: a TaxTree structure.
        lineage_file: name of the lineage file
        nodes: a counter for TaxIds
        collapse: This bool controls the collapse of taxid 131567
            (cellular organisms) and is True by default

    Returns: A string with the output messages

    """
    log, taxids_dic = tree.get_lineage(parents, iter(nodes))
    output: io.StringIO = io.StringIO(newline='')
    output.write(log)
    if collapse:  # Collapse taxid 131567 (cellular organisms) if desired
        for tid in taxids_dic:
            if len(taxids_dic[tid]) > 2:  # Not collapse for unclassified
                try:
                    taxids_dic[tid].remove(CELLULAR_ORGANISMS)
                except ValueError:
                    pass
    lineage_dic = {tax_id: [names[tid] for tid in taxids_dic[tax_id]]
                   for tax_id in taxids_dic}
    output.write(f'  \033[90mSaving lineage file {lineage_file} with '
                 f'{len(nodes)} nodes...\033[0m')
    with open(lineage_file, 'w', newline='') as tsv_handle:
        tsvwriter = csv.writer(tsv_handle, dialect='krona')
        tsvwriter.writerow(["#taxID", "Lineage"])
        for tid in nodes:
            counts: str = str(nodes[tid])  # nodes[tid] is a int
            row: List[Union[TaxId, str]] = [counts, ]
            row.extend(lineage_dic[tid])
            tsvwriter.writerow(row)
    output.write('\033[92m OK! \033[0m\n')
    return output.getvalue()


def krona_from_text(samples: List[Sample],
                    outputs: Dict[Rank, Filename],
                    htmlfile: Filename = Filename('Output' + HTML_SUFFIX),
                    ):
    """Generate the Krona html file calling ktImportText.

    Superseded by krona.krona_from_xml().

    """
    subprc = ["ktImportText"]
    subprc.extend(samples)
    try:
        subprc.extend([outputs[level][i]
                       for level in List(Rank.selected_ranks)
                       for i in range(len(outputs[level]))])
    except KeyError:
        pass
    subprc.extend(["-o", htmlfile])
    try:
        subprocess.run(subprc, check=True)
    except subprocess.CalledProcessError:
        print('\n\033[91mERROR!\033[0m ktImportText: ' +
              'returned a non-zero exit status (Krona plot built failed)')
