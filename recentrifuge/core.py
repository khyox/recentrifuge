"""
Recentrifuge core functions.

"""
import collections as col
import csv
import io
import statistics
import subprocess
import sys
from typing import List, Set, Tuple, Union, Dict, Counter, Optional

from recentrifuge.config import Filename, Sample, Id, Parents, Score, Scores
from recentrifuge.config import HTML_SUFFIX, CELLULAR_ORGANISMS
from recentrifuge.config import STR_CONTROL, STR_EXCLUSIVE, STR_SHARED
from recentrifuge.config import STR_SUMMARY, STR_CONTROL_SHARED
from recentrifuge.config import UnionCounter, UnionScores
from recentrifuge.config import gray, red, yellow, blue, magenta, cyan, green
from recentrifuge.ontology import Ontology
from recentrifuge.params import EPS, ROBUST_MIN_SAMPLES, ROBUST_XOVER_OUTLIER
from recentrifuge.params import MILD_CONTM_MIN_RELFREQ
from recentrifuge.params import ROBUST_XOVER_ORD_MAG, SEVR_CONTM_MIN_RELFREQ
from recentrifuge.rank import Rank, TaxLevels
from recentrifuge.shared_counter import SharedCounter
from recentrifuge.trees import TaxTree, SampleDataById


def process_rank(*args,
                 **kwargs
                 ) -> Tuple[List[Sample],
                            Dict[Sample, UnionCounter],
                            Dict[Sample, Counter[Id]],
                            Dict[Sample, UnionScores]]:
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """

    # Recover input and parameters
    rank: Rank = args[0]
    controls: int = kwargs['controls']
    mintaxas: Dict[Sample, int] = kwargs['mintaxas']
    ontology: Ontology = kwargs['ontology']
    including = ontology.including
    excluding = ontology.excluding
    taxids: Dict[Sample, TaxLevels] = kwargs['taxids']
    counts: Dict[Sample, UnionCounter] = kwargs['counts']
    accs: Dict[Sample, Counter[Id]] = kwargs['accs']
    scores: Dict[Sample, UnionScores] = kwargs['scores']
    raws: List[Sample] = kwargs['raw_samples']
    output: io.StringIO = io.StringIO(newline='')

    def vwrite(*args) -> None:
        """Print only if verbose/debug mode is enabled"""
        if kwargs['debug']:
            output.write(' '.join(str(item) for item in args))

    def fltlst2str(lst: List[float]) -> str:
        """Convert a list of floats into a nice string"""
        return '[' + gray((', '.join(f'{elm:.1g}' for elm in lst))) + ']'

    def blst2str(lst: List[bool]) -> str:
        """Convert a list of booleans into a nice string"""
        return ('[' + (', '.join(magenta('T') if elm else 'F' for elm in lst))
                + ']')

    def get_shared_mintaxa() -> int:
        """Give a value of mintaxa for shared derived samples

        This value is currently the minimum of the mintaxa of all the
         (non control) raw samples.
        """
        return min([mintaxas[smpl] for smpl in raws[controls:]])

    # Declare/define variables
    samples: List[Sample] = []
    # pylint: disable = unused-variable
    shared_counts: SharedCounter = SharedCounter()
    shared_score: SharedCounter = SharedCounter()
    shared_ctrl_counts: SharedCounter = SharedCounter()
    shared_ctrl_score: SharedCounter = SharedCounter()
    # pylint: enable = unused-variable

    output.write(f'\033[90mAnalysis for taxonomic rank "'
                 f'\033[95m{rank.name.lower()}\033[90m":\033[0m\n')

    def cross_analysis(iteration, raw):
        """Cross analysis: exclusive and part of shared&ctrl"""
        nonlocal shared_counts, shared_score
        nonlocal shared_ctrl_counts, shared_ctrl_score

        def partial_shared_update(i):
            """Perform shared and shared-control taxa partial evaluations"""
            nonlocal shared_counts, shared_score
            nonlocal shared_ctrl_counts, shared_ctrl_score
            if i == 0:  # 1st iteration: Initialize shared abundance and score
                shared_counts.update(sub_shared_counts)
                shared_score.update(sub_shared_score)
            elif i < controls:  # Just update shared abundance and score
                shared_counts &= sub_shared_counts
                shared_score &= sub_shared_score
            elif i == controls:  # Initialize shared-control counters
                shared_counts &= sub_shared_counts
                shared_score &= sub_shared_score
                shared_ctrl_counts.update(sub_shared_counts)
                shared_ctrl_score.update(sub_shared_score)
            elif controls:  # Both: Accumulate shared abundance and score
                shared_counts &= sub_shared_counts
                shared_score &= sub_shared_score
                shared_ctrl_counts &= sub_shared_counts
                shared_ctrl_score &= sub_shared_score
            else:  # Both: Accumulate shared abundance and score (no controls)
                shared_counts &= sub_shared_counts
                shared_score &= sub_shared_score

        exclude: Set[Id] = set()
        # Get taxids at this rank that are present in the other samples
        for sample in (smpl for smpl in raws if smpl != raw):
            exclude.update(taxids[sample][rank])
        exclude.update(excluding)  # Add explicit excluding taxa if any
        output.write(f'  \033[90mExclusive: From \033[0m{raw}\033[90m '
                     f'excluding {len(exclude)} taxa. '
                     f'Generating sample...\033[0m')

        exclude_tree = TaxTree()
        exclude_out = SampleDataById(['counts', 'scores', 'accs'])
        exclude_tree.allin1(ontology=ontology,
                            counts=counts[raw],
                            scores=scores[raw],
                            min_taxa=mintaxas[raw],
                            min_rank=rank,
                            just_min_rank=True,
                            include=including,
                            exclude=exclude,
                            out=exclude_out)
        exclude_out.purge_counters()
        if exclude_out.counts:  # Avoid adding empty samples
            sample = Sample(f'{raw}_{STR_EXCLUSIVE}_{rank.name.lower()}')
            samples.append(sample)
            counts[sample] = exclude_out.get_counts()
            accs[sample] = exclude_out.get_accs()
            scores[sample] = exclude_out.get_scores()
            output.write('\033[92m OK! \033[0m\n')
        else:
            output.write('\033[93m VOID \033[0m\n')

        # Get partial abundance and score for the shared analysis
        sub_shared_tree = TaxTree()
        sub_shared_out = SampleDataById(['shared', 'accs'])
        sub_shared_tree.allin1(ontology=ontology,
                               counts=counts[raw],
                               scores=scores[raw],
                               min_taxa=mintaxas[raw],
                               min_rank=rank,
                               just_min_rank=True,
                               include=including,
                               exclude=excluding,
                               out=sub_shared_out)
        sub_shared_out.purge_counters()
        # Scale scores by abundance
        sub_shared_counts: SharedCounter = sub_shared_out.get_shared_counts()
        sub_shared_score: SharedCounter = sub_shared_out.get_shared_scores()
        sub_shared_score *= sub_shared_counts
        partial_shared_update(iteration)

    def shared_analysis():
        """Perform last steps of shared taxa analysis"""
        shared_tree: TaxTree = TaxTree()
        shared_out: SampleDataById = SampleDataById(['shared', 'accs'])
        shared_tree.allin1(ontology=ontology,
                           counts=shared_counts,
                           scores=shared_score,
                           min_taxa=get_shared_mintaxa(),
                           include=including,
                           exclude=excluding,
                           out=shared_out)
        shared_out.purge_counters()
        out_counts: SharedCounter = shared_out.get_shared_counts()
        output.write(gray(f'  Shared: Including {len(out_counts)}'
                          ' shared taxa. Generating sample... '))
        if out_counts:
            sample = Sample(f'{STR_SHARED}_{rank.name.lower()}')
            samples.append(sample)
            counts[Sample(sample)] = out_counts
            accs[Sample(sample)] = shared_out.get_accs()
            scores[sample] = shared_out.get_shared_scores()
            output.write(green('OK!\n'))
        else:
            output.write(yellow('VOID\n'))

    def control_analysis():
        """Perform last steps of control and shared controls analysis"""
        nonlocal shared_ctrl_counts, shared_ctrl_score

        def robust_contamination_removal():
            """Implement robust contamination removal algorithm."""
            nonlocal exclude_sets, shared_crossover

            def compute_qn(data: List[float], dist: str = "Gauss") -> float:
                """Compute Qn robust estimator of scale (Rousseeuw, 1993)"""
                c_d: float  # Select d parameter depending on the distribution
                if dist == "Gauss":
                    c_d = 2.2219
                elif dist == "Cauchy":  # Heavy-tailed distribution
                    c_d = 1.2071
                elif dist == "NegExp":  # Negative exponential (asymetric)
                    c_d = 3.4760
                else:
                    raise Exception(red('\nERROR! ') + 'Unknown distribution')
                num: int = len(data)
                sort_data = sorted(data)
                pairwisedifs: List[float] = []
                for (i, x_val) in enumerate(sort_data):
                    for y_val in sort_data[i + 1:]:
                        pairwisedifs.append(abs(x_val - y_val))
                k: int = int(num * (num / 2 + 1) / 4)
                return c_d * sorted(pairwisedifs)[k - 1]

            exclude_sets = {smpl: set() for smpl in raws[controls:]}
            vwrite(gray('Robust contamination removal: '
                        'Searching for contaminants...\n'))
            for tid in exclude_candidates:
                relfreq_ctrl: List[float] = [accs[ctrl][tid]
                                             / accs[ctrl][ontology.ROOT]
                                             for ctrl in raws[:controls]]
                relfreq_smpl: List[float] = [accs[smpl][tid]
                                             / accs[smpl][ontology.ROOT]
                                             for smpl in raws[controls:]]
                relfreq: List[float] = relfreq_ctrl + relfreq_smpl
                crossover: List[bool]  # Crossover source (yes/no)
                # Just-controls contamination check
                if all([rf < EPS for rf in relfreq_smpl]):
                    vwrite(cyan('just-ctrl:\t'), tid, ontology.get_name(tid),
                           gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                           fltlst2str(relfreq_smpl), '\n')
                    continue  # Go for next candidate
                # Critical contamination check
                if all([rf > SEVR_CONTM_MIN_RELFREQ for rf in relfreq_ctrl]):
                    vwrite(red('critical:\t'), tid, ontology.get_name(tid),
                           gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                           fltlst2str(relfreq_smpl), '\n')
                    for exclude_set in exclude_sets.values():
                        exclude_set.add(tid)
                    continue  # Go for next candidate
                # Severe contamination check
                if any([rf > SEVR_CONTM_MIN_RELFREQ for rf in relfreq_ctrl]):
                    vwrite(yellow('severe: \t'), tid, ontology.get_name(tid),
                           gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                           fltlst2str(relfreq_smpl), '\n')
                    for exclude_set in exclude_sets.values():
                        exclude_set.add(tid)
                    continue  # Go for next candidate
                # Mild contamination check
                if all([rf > MILD_CONTM_MIN_RELFREQ for rf in relfreq_ctrl]):
                    vwrite(blue('mild cont:\t'), tid, ontology.get_name(tid),
                           gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                           fltlst2str(relfreq_smpl), '\n')
                    for exclude_set in exclude_sets.values():
                        exclude_set.add(tid)
                    continue  # Go for next candidate
                # Calculate median and MAD median but including controls
                mdn: float = statistics.median(relfreq)
                # mad:float=statistics.mean([abs(mdn - rf) for rf in relfreq])
                q_n: float = compute_qn(relfreq, dist="NegExp")
                # Calculate crossover in samples
                outlier_lim: float = mdn + ROBUST_XOVER_OUTLIER * q_n
                ordomag_lim: float = max(
                    relfreq_ctrl) * 10 ** ROBUST_XOVER_ORD_MAG
                crossover = [rf > outlier_lim and rf > ordomag_lim
                             for rf in relfreq[controls:]]
                # Crossover contamination check
                if any(crossover):
                    vwrite(magenta('crossover:\t'), tid,
                           ontology.get_name(tid), green(
                            f'lims: [{outlier_lim:.1g}]' + (
                                '<' if outlier_lim < ordomag_lim else '>') +
                            f'[{ordomag_lim:.1g}]'),
                           gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                           fltlst2str(relfreq_smpl),
                           gray('crossover:'), blst2str(crossover), '\n')
                    # Exclude just for contaminated samples (not the source)
                    vwrite(magenta('\t->'), gray(f'Include {tid} just in:'))
                    for i in range(len(raws[controls:])):
                        if not crossover[i]:
                            exclude_sets[raws[i + controls]].add(tid)
                        else:
                            vwrite(f' {raws[i + controls]}')
                    if all(crossover):  # Shared taxon contaminating control(s)
                        vwrite(' (', yellow('Shared crossover taxon!'), ')')
                        shared_crossover.add(tid)
                    vwrite('\n')
                    continue
                # Other contamination: remove from all samples
                vwrite(gray('other cont:\t'), tid, ontology.get_name(tid),
                       green(f'lims: [{outlier_lim:.1g}]' + (
                           '<' if outlier_lim < ordomag_lim else '>') +
                             f'[{ordomag_lim:.1g}]'),
                       gray('relfreq:'), fltlst2str(relfreq_ctrl) +
                       fltlst2str(relfreq_smpl), '\n')
                for exclude_set in exclude_sets.values():
                    exclude_set.add(tid)

        # Get taxids at this rank that are present in the control samples
        exclude_candidates: Set[Id] = set()
        for i in range(controls):
            exclude_candidates.update(taxids[raws[i]][rank])
        exclude_sets: Dict[Sample, Set[Id]]
        shared_crossover: Set[Id] = set()  # Shared taxa contaminating controls
        if controls and (len(raws) - controls >= ROBUST_MIN_SAMPLES):
            robust_contamination_removal()
        else:  # If this case, just apply strict control
            exclude_sets = {file: exclude_candidates
                            for file in raws[controls::]}
        # Add explicit excluding taxa (if any) to exclude sets
        for exclude_set in exclude_sets.values():
            exclude_set.update(excluding)
        exclude_candidates.update(excluding)
        # Process each sample excluding control taxa
        for raw in raws[controls:]:
            output.write(gray('  Ctrl: From') + f' {raw} ' +
                         gray(f'excluding {len(exclude_sets[raw])} ctrl taxa. '
                              f'Generating sample... '))
            ctrl_tree = TaxTree()
            ctrl_out = SampleDataById(['counts', 'scores', 'accs'])
            ctrl_tree.allin1(ontology=ontology,
                             counts=counts[raw],
                             scores=scores[raw],
                             min_taxa=mintaxas[raw],
                             min_rank=rank,
                             just_min_rank=True,
                             include=including,
                             exclude=exclude_sets[raw],
                             out=ctrl_out)
            ctrl_out.purge_counters()
            if ctrl_out.counts:  # Avoid adding empty samples
                sample = Sample(f'{raw}_{STR_CONTROL}_{rank.name.lower()}')
                samples.append(sample)
                counts[sample] = ctrl_out.get_counts()
                accs[sample] = ctrl_out.get_accs()
                scores[sample] = ctrl_out.get_scores()
                output.write(green('OK!\n'))
            else:
                output.write(yellow('VOID\n'))

        def shared_ctrl_analysis():
            """Perform last steps of shared taxa analysis"""
            shared_ctrl_tree: TaxTree = TaxTree()
            shared_ctrl_out: SampleDataById = SampleDataById(
                ['shared', 'accs'])
            shared_ctrl_tree.allin1(ontology=ontology,
                                    counts=shared_ctrl_counts,
                                    scores=shared_ctrl_score,
                                    min_taxa=get_shared_mintaxa(),
                                    include=including,
                                    exclude=(exclude_candidates
                                             - shared_crossover),
                                    out=shared_ctrl_out)
            shared_ctrl_out.purge_counters()
            out_counts: SharedCounter = shared_ctrl_out.get_shared_counts()
            output.write(gray(f'  Ctrl-shared: Including {len(out_counts)}'
                              ' shared taxa. Generating sample... '))
            if out_counts:
                sample = Sample(f'{STR_CONTROL_SHARED}_{rank.name.lower()}')
                samples.append(sample)
                counts[Sample(sample)] = out_counts
                accs[Sample(sample)] = shared_ctrl_out.get_accs()
                scores[sample] = shared_ctrl_out.get_shared_scores()
                output.write(green('OK!\n'))
            else:
                output.write(yellow('VOID\n'))

        # Shared-control taxa final analysis
        if shared_ctrl_counts:
            # Normalize scaled scores by total abundance
            shared_ctrl_score /= (+shared_ctrl_counts)
            # Get averaged abundance by number of samples minus ctrl samples
            shared_ctrl_counts //= (len(raws) - controls)
            shared_ctrl_analysis()
        else:
            output.write(gray('  Ctrl-shared: No taxa! ') +
                         yellow('VOID') + gray(' sample.\n'))

    # Cross analysis iterating by output: exclusive and part of shared&ctrl
    for num_file, raw_sample_name in enumerate(raws):
        cross_analysis(num_file, raw_sample_name)

    # Shared taxa final analysis
    shared_counts = +shared_counts  # remove counts <= 0
    if shared_counts:
        # Normalize scaled scores by total abundance (after eliminating zeros)
        shared_score /= (+shared_counts)
        # Get averaged abundance by number of samples
        shared_counts //= len(raws)
        shared_analysis()
    else:
        output.write(gray('  Shared: No shared taxa! ') +
                     yellow('VOID') + gray(' sample.\n'))

    # Control sample subtraction
    if controls:
        control_analysis()

    # Print output and return
    print(output.getvalue())
    sys.stdout.flush()
    return samples, counts, accs, scores


def summarize_analysis(*args,
                       **kwargs
                       ) -> Tuple[Optional[Sample],
                                  Counter[Id],
                                  Counter[Id],
                                  Scores]:
    """
    Summarize for a cross-analysis (to be usually called in parallel!).
    """
    # Recover input and parameters
    analysis: str = args[0]
    ontology: Ontology = kwargs['ontology']
    # TODO: Delete the following comment lines in a future release
    # including = ontology.including   # See comment below for the reason
    # excluding = ontology.excluding   # in/excluding are not used anymore
    counts: Dict[Sample, Counter[Id]] = kwargs['counts']
    scores: Dict[Sample, Dict[Id, Score]] = kwargs['scores']
    samples: List[Sample] = kwargs['samples']
    output: io.StringIO = io.StringIO(newline='')

    # Declare/define variables
    summary_counts: Counter[Id] = col.Counter()
    summary_acc: Counter[Id] = col.Counter()
    summary_score: Scores = Scores({})
    summary: Optional[Sample] = None

    output.write(gray('Summary for ') + analysis + gray('... '))

    target_samples: List[Sample] = [smpl for smpl in samples
                                    if smpl.startswith(analysis)]
    assert len(target_samples) >= 1, \
        red('ERROR! ') + analysis + gray(' has no samples to summarize!')
    for smpl in target_samples:
        summary_counts += counts[smpl]
        summary_score.update(scores[smpl])

    tree = TaxTree()
    tree.grow(ontology=ontology,
              counts=summary_counts,
              scores=summary_score)
    tree.subtract()
    tree.shape()
    summary_counts.clear()
    summary_score.clear()
    # Avoid including/excluding here as get_taxa is not as 'clever' as allin1
    #  and taxa are already included/excluded in the derived samples
    tree.get_taxa(counts=summary_counts,
                  accs=summary_acc,
                  scores=summary_score)
    summary_counts = +summary_counts  # remove counts <= 0
    if summary_counts:  # Avoid returning empty sample (summary would be None)
        summary = Sample(f'{analysis}_{STR_SUMMARY}')
        output.write(gray('(') + cyan(f'{len(target_samples)}') +
                     gray(' samples)') + green(' OK!\n'))
    else:
        output.write(yellow(' VOID\n'))
    # Print output and return
    print(output.getvalue(), end='')
    sys.stdout.flush()
    return summary, summary_counts, summary_acc, summary_score


def write_lineage(ontology: Ontology,
                  parents: Parents,
                  names: Dict[Id, str],
                  tree: TaxTree,
                  lineage_file: str,
                  nodes: Counter[Id],
                  collapse: bool = True,
                  ) -> str:
    """
    Writes a lineage file understandable by Krona.

    Args:
        ontology: any Ontology object.
        parents: dictionary of parents for every Id.
        names: dictionary of names for every Id.
        tree: a TaxTree structure.
        lineage_file: name of the lineage file
        nodes: a counter for TaxIds
        collapse: This bool controls the collapse of taxid 131567
            (cellular organisms) and is True by default

    Returns: A string with the output messages

    """
    log, taxids_dic = tree.get_lineage(ontology, parents, iter(nodes))
    output: io.StringIO = io.StringIO(newline='')
    output.write(log)
    # TODO: Generalize this for non-NCBI taxonomies (collapse is specific)
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
            row: List[Union[Id, str]] = [counts, ]
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
                       for level in list(Rank.selected_ranks)
                       for i in range(len(outputs[level]))])
    except KeyError:
        pass
    subprc.extend(["-o", htmlfile])
    try:
        subprocess.run(subprc, check=True)
    except subprocess.CalledProcessError:
        print('\n\033[91mERROR!\033[0m ktImportText: ' +
              'returned a non-zero exit status (Krona plot built failed)')
