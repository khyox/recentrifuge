#!/usr/bin/env python3
"""
Post-process Centrifuge/LMAT output.
"""
# pylint: disable=no-name-in-module, not-an-iterable
import argparse
import collections as col
import multiprocessing as mp
import os
import platform
import sys
import time
from typing import Counter, List, Dict, Set, Callable, Tuple

from recentrifuge.centrifuge import process_report, process_output
from recentrifuge.centrifuge import select_centrifuge_inputs
from recentrifuge.config import Filename, Sample, TaxId, Score, Scoring, Excel
from recentrifuge.config import HTML_SUFFIX, DEFMINTAXA, TAXDUMP_PATH
from recentrifuge.config import NODES_FILE, NAMES_FILE, PLASMID_FILE
from recentrifuge.config import STR_CONTROL, STR_EXCLUSIVE, STR_SHARED
from recentrifuge.config import STR_CONTROL_SHARED, Err
from recentrifuge.config import gray, red, green, yellow, blue, magenta
from recentrifuge.core import process_rank, summarize_analysis
from recentrifuge.krona import COUNT, UNASSIGNED, SCORE
from recentrifuge.krona import KronaTree
from recentrifuge.lmat import select_lmat_inputs
from recentrifuge.rank import Rank, TaxLevels
from recentrifuge.taxonomy import Taxonomy
from recentrifuge.trees import TaxTree, MultiTree, SampleDataByTaxId

# optional package pandas (to generate Excel output)
_USE_PANDAS = True
try:
    import pandas as pd
except ImportError:
    pd = None
    _USE_PANDAS = False

__version__ = '0.17.3_rc2'
__author__ = 'Jose Manuel Marti'
__date__ = 'Feb 2018'


def _debug_dummy_plot(taxonomy: Taxonomy,
                      htmlfile: Filename,
                      scoring: Scoring = Scoring.SHEL,
                      ):
    """

    Generate dummy Krona plot via Krona 2.0 XML spec and exit

    """
    print(gray(f'Generating dummy Krona plot {htmlfile}...'), end='')
    sys.stdout.flush()
    samples: List[Sample] = [Sample('SINGLE'), ]
    krona: KronaTree = KronaTree(samples,
                                 min_score=Score(35),
                                 max_score=Score(100),
                                 scoring=scoring,
                                 )
    polytree: MultiTree = MultiTree(samples=samples)
    polytree.grow(taxonomy=taxonomy)
    polytree.toxml(taxonomy=taxonomy, krona=krona)
    krona.tohtml(htmlfile, pretty=True)
    print(green('OK!'))


def main():
    """Main entry point to Recentrifuge."""

    def configure_parser():
        """Argument Parser Configuration"""
        parser = argparse.ArgumentParser(
            description='Post-process Centrifuge/LMAT output',
            epilog=f'%(prog)s  - {__author__} - {__date__}',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            '-n', '--nodespath',
            action='store',
            metavar='PATH',
            default=TAXDUMP_PATH,
            help=('path for the nodes information files '
                  '(nodes.dmp and names.dmp from NCBI)')
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            version=f'%(prog)s release {__version__} ({__date__})'
        )
        parser_filein = parser.add_mutually_exclusive_group(required=True)
        parser_filein.add_argument(
            '-f', '--file',
            action='append',
            metavar='FILE',
            type=Filename,
            help=('Centrifuge output files. If a single directory is entered, '
                  'every .out file inside will be taken as a different sample.'
                  ' Multiple -f is available to include several samples.')
        )
        parser_filein.add_argument(
            '-r', '--report',
            action='append',
            metavar='FILE',
            type=Filename,
            help=('Centrifuge/Kraken report files '
                  '(multiple -r is available to include several samples)')
        )
        parser_filein.add_argument(
            '-l', '--lmat',
            action='append',
            metavar='FILE',
            type=Filename,
            default=None,
            help=('LMAT output dir or file prefix. If just "." is entered, '
                  'every subdirectory under the current directory will be '
                  'taken as a sample and scanned looking for LMAT output files'
                  '. Multiple -l is available to include several samples.')
        )
        parser_cross = parser.add_mutually_exclusive_group(required=False)
        parser_cross.add_argument(
            '-a', '--avoidcross',
            action='store_true',
            help='avoid cross analysis'
        )
        parser_cross.add_argument(
            '-c', '--controls',
            action='store',
            metavar='CONTROLS_NUMBER',
            type=int,
            default=0,
            help=('this number of first samples will be treated as negative '
                  'controls; default is no controls')
        )
        parser_output = parser.add_argument_group(
            'output', 'Related to the output files')
        parser_output.add_argument(
            '-e', '--excel',
            action='store',
            metavar='OUTPUT_TYPE',
            choices=[str(excel) for excel in Excel],
            default=str(Excel(0)),
            help=(f'type of scoring to be applied, and can be one of '
                  f'{[str(excel) for excel in Excel]}')
        )
        parser_output.add_argument(
            '-o', '--outhtml',
            action='store',
            metavar='FILE',
            type=Filename,
            help='HTML output file (if not given the filename will be '
                 'inferred from input files)'
        )
        parser_mode = parser.add_argument_group('mode',
                                                'Specific modes of running')
        parser_mode.add_argument(
            '--dummy',  # hidden flag: just generate a dummy plot for JS debug
            action='store_true',
            help=argparse.SUPPRESS
        )
        parser_mode.add_argument(
            '-g', '--debug',
            action='store_true',
            help='increase output verbosity and perform additional checks'
        )
        parser_mode.add_argument(
            '--sequential',
            action='store_true',
            help='deactivate parallel processing'
        )
        parser_tuning = parser.add_argument_group(
            'tuning',
            'Fine tuning of algorithm parameters')
        parser_tuning.add_argument(
            '-i', '--include',
            action='append',
            metavar='TAXID',
            type=TaxId,
            default=[],
            help=('NCBI taxid code to include a taxon and all underneath '
                  '(multiple -i is available to include several taxid); '
                  'by default all the taxa is considered for inclusion')
        )
        parser_tuning.add_argument(
            '-k', '--nokollapse',
            action='store_true',
            help='show the "cellular organisms" taxon'
        )
        parser_tuning.add_argument(
            '-m', '--mintaxa',
            action='store',
            metavar='INT',
            type=int,
            default=DEFMINTAXA,
            help='minimum taxa to avoid collapsing one level to the parent one'
        )
        parser_tuning.add_argument(
            '-s', '--scoring',
            action='store',
            metavar='SCORING',
            choices=[str(each_score) for each_score in Scoring],
            default=str(Scoring(0)),
            help=(f'type of scoring to be applied, and can be one of '
                  f'{[str(scoring) for scoring in Scoring]}')
        )
        parser_tuning.add_argument(
            '-t', '--takeoutroot',
            action='store_true',
            help='remove counts directly assigned to the "root" level'
        )
        parser_tuning.add_argument(
            '-u', '--summary',
            action='store',
            metavar='OPTION',
            choices=['add', 'only', 'avoid'],
            default='add',
            help=('select to "add" summary samples to other samples, or to '
                  '"only" show summary samples or to "avoid" summaries at all')
        )
        parser_tuning.add_argument(
            '-w', '--ctrlmintaxa',
            action='store',
            metavar='INT',
            type=int,
            default=None,
            help='minimum taxa to avoid collapsing one level to the parent one'
                 ' in control samples; it defaults to "mintaxa"'
        )
        parser_tuning.add_argument(
            '-x', '--exclude',
            action='append',
            metavar='TAXID',
            type=TaxId,
            default=[],
            help=('NCBI taxid code to exclude a taxon and all underneath '
                  '(multiple -x is available to exclude several taxid)')
        )
        parser_tuning.add_argument(
            '-y', '--minscore',
            action='store',
            metavar='NUMBER',
            type=lambda txt: Score(float(txt)),
            default=None,
            help=('minimum score/confidence of the classification of a read '
                  'to pass the quality filter; all pass by default')
        )
        parser_tuning.add_argument(
            '-z', '--ctrlminscore',
            action='store',
            metavar='NUMBER',
            type=lambda txt: Score(float(txt)),
            default=None,
            help=('minimum score/confidence of the classification of a read '
                  'in control samples to pass the quality filter; if defaults '
                  'to "minscore"')
        )
        return parser

    def check_debug():
        """Check debugging mode"""
        if args.debug:
            print(gray('INFO: Debugging mode activated\n'))

    def select_inputs():
        """Choose right input and output files"""
        nonlocal process, scoring, input_files, plasmidfile

        if outputs and len(outputs) == 1 and os.path.isdir(outputs[0]):
            select_centrifuge_inputs(outputs)
        if lmats:
            plasmidfile = Filename(os.path.join(args.nodespath, PLASMID_FILE))
            select_lmat_inputs(lmats)

        # Select method and arguments depending on type of files to analyze
        if lmats:
            process = process_output
            input_files = lmats
            scoring = Scoring.LMAT
        elif reports:
            process = process_report
            input_files = reports
        else:
            process = process_output
            input_files = outputs

    def check_controls():
        """Check and info about the control samples"""
        if args.controls:
            if args.controls > len(input_files):
                print(red(' ERROR!'), gray('More controls than samples'))
                exit(1)
            print(gray('Control(s) sample(s) for subtractions:'))
            for i in range(args.controls):
                print(blue(f'\t{input_files[i]}'))

    def select_html_file():
        """HTML filename selection"""
        nonlocal htmlfile
        if lmats:  # Select case for dir name or filename prefix
            if os.path.isdir(lmats[0]):  # Dir name
                dirname = os.path.dirname(os.path.normpath(lmats[0]))
                if not dirname or dirname == '.':
                    basename = 'output'
                else:
                    basename = os.path.basename(dirname)
            else:  # Explicit path and file name prefix is provided
                dirname, basename = os.path.split(lmats[0])
            htmlfile = Filename(os.path.join(dirname, basename + HTML_SUFFIX))
        elif reports:
            htmlfile = Filename(reports[0].split('_mhl')[0] + HTML_SUFFIX)
        else:
            htmlfile = Filename(outputs[0].split('_mhl')[0] + HTML_SUFFIX)

    def read_samples():
        """Read samples"""
        print(gray('\nPlease, wait, processing files in parallel...\n'))
        # Enable parallelization with 'spawn' under known platforms
        if platform.system() and not args.sequential:  # Only for known systems
            mpctx = mp.get_context('fork')
            with mpctx.Pool(processes=min(os.cpu_count(),
                                          len(input_files))) as pool:
                async_results = [pool.apply_async(
                    process,
                    args=[input_files[num],  # file name
                          True if num < args.controls else False],  # is ctrl?
                    kwds=kwargs
                ) for num in range(len(input_files))]
                for file, (sample, tree, out, err) in zip(
                        input_files, [r.get() for r in async_results]):
                    if err is Err.NO_ERROR:
                        samples.append(sample)
                        trees[sample] = tree
                        taxids[sample] = out.get_taxlevels()
                        counts[sample] = out.counts
                        accs[sample] = out.accs
                        scores[sample] = out.scores
                    elif err is Err.VOID_CTRL:
                        print('There were void controls.', red('Aborting!'))
                        exit(1)
        else:  # sequential processing of each sample
            for num, file in enumerate(input_files):
                (sample, tree, out, err) = process(
                    file, True if num < args.controls else False, **kwargs)
                if err is Err.NO_ERROR:
                    samples.append(sample)
                    trees[sample] = tree
                    taxids[sample] = out.get_taxlevels()
                    counts[sample] = out.counts
                    accs[sample] = out.accs
                    scores[sample] = out.scores
                elif err is Err.VOID_CTRL:
                    print('There were void controls.', red('Aborting!'))
                    exit(1)
        raw_samples.extend(samples)  # Store raw sample names

    def analyze_samples():
        """Cross analysis of samples in parallel by taxlevel"""
        print(gray('Please, wait. Performing cross analysis in parallel...\n'))
        # Update kwargs with more parameters for the followings func calls
        kwargs.update({'taxids': taxids, 'counts': counts, 'scores': scores,
                       'accs': accs, 'raw_samples': raw_samples})
        if platform.system() and not args.sequential:  # Only for known systems
            mpctx = mp.get_context('fork')  # Important for OSX&Win
            with mpctx.Pool(processes=min(os.cpu_count(), len(
                    Rank.selected_ranks))) as pool:
                async_results = [pool.apply_async(
                    process_rank,
                    args=[level],
                    kwds=kwargs
                ) for level in Rank.selected_ranks]
                for level, (smpls, abunds, accumulators, score) in zip(
                        Rank.selected_ranks,
                        [r.get() for r in async_results]):
                    samples.extend(smpls)
                    counts.update(abunds)
                    accs.update(accumulators)
                    scores.update(score)
        else:  # sequential processing of each selected rank
            for level in Rank.selected_ranks:
                (smpls, abunds,
                 accumulators, score) = process_rank(level, **kwargs)
                samples.extend(smpls)
                counts.update(abunds)
                accs.update(accumulators)
                scores.update(score)

    def summarize_samples():
        """Summary of samples in parallel by type of cross-analysis"""
        # timing initialization
        summ_start_time: float = time.perf_counter()
        print(gray('Please, wait. Generating summaries in parallel...'))
        # Update kwargs with more parameters for the followings func calls
        kwargs.update({'samples': samples})
        # Get list of set of samples to summarize (note pylint bug #776)
        # pylint: disable=unsubscriptable-object
        target_analysis: col.OrderedDict[str, None] = col.OrderedDict({
            f'{raw}_{study}': None for study in [STR_EXCLUSIVE, STR_CONTROL]
            for raw in raw_samples
            for smpl in samples if smpl.startswith(f'{raw}_{study}')
        })
        # pylint: enable=unsubscriptable-object
        # Add shared and control_shared analysis if they exist (are not void)
        for study in [STR_SHARED, STR_CONTROL_SHARED]:
            for smpl in samples:
                if smpl.startswith(study):
                    target_analysis[study] = None
                    break

        if platform.system() and not args.sequential:  # Only for known systems
            mpctx = mp.get_context('fork')
            with mpctx.Pool(processes=min(os.cpu_count(),
                                          len(input_files))) as pool:
                async_results = [pool.apply_async(
                    summarize_analysis,
                    args=[analysis],
                    kwds=kwargs
                ) for analysis in target_analysis]
                for analysis, (summary, abund, acc, score) in zip(
                        target_analysis, [r.get() for r in async_results]):
                    if summary:  # Avoid adding empty samples
                        summaries.append(summary)
                        counts[summary] = abund
                        accs[summary] = acc
                        scores[summary] = score
        else:  # sequential processing of each selected rank
            for analysis in target_analysis:
                (summary, abund,
                 acc, score) = summarize_analysis(analysis, **kwargs)
                if summary:  # Avoid adding empty samples
                    summaries.append(summary)
                    counts[summary] = abund
                    accs[summary] = acc
                    scores[summary] = score
        # Timing results
        print(gray('Summary elapsed time:'),
              f'{time.perf_counter() - summ_start_time:.3g}', gray('sec'))

    def generate_krona():
        """Generate Krona plot with all the results via Krona 2.0 XML spec"""

        print(gray('\nBuilding the taxonomy multiple tree... '), end='')
        sys.stdout.flush()
        krona: KronaTree = KronaTree(samples,
                                     num_raw_samples=len(raw_samples),
                                     min_score=Score(
                                         min([min(scores[sample].values())
                                              for sample in samples
                                              if len(scores[sample])])),
                                     max_score=Score(
                                         max([max(scores[sample].values())
                                              for sample in samples
                                              if len(scores[sample])])),
                                     scoring=scoring,
                                     )
        polytree.grow(taxonomy=ncbi,
                      abundances=counts,
                      accs=accs,
                      scores=scores)
        print(green('OK!'))
        print(gray('Generating final plot (') + magenta(htmlfile) +
              gray(')... '), end='')
        sys.stdout.flush()
        polytree.toxml(taxonomy=ncbi, krona=krona)
        krona.tohtml(htmlfile, pretty=False)
        print(green('OK!'))

    def generate_excel():
        """Generate Excel with results via pandas DataFrame"""

        xlsx_name: Filename = Filename(htmlfile.split('.html')[0] + '.xlsx')
        print(gray(f'Generating Excel {str(excel).lower()} summary (') +
              magenta(xlsx_name) + gray(')... '), end='')
        sys.stdout.flush()
        xlsxwriter = pd.ExcelWriter(xlsx_name)
        list_rows: List = []
        data_frame: pd.DataFrame
        if excel is Excel.FULL:
            polytree.to_items(taxonomy=ncbi, items=list_rows)
            # Generate the pandas DataFrame from items and export to Excel
            iterable_1 = [samples, [COUNT, UNASSIGNED, SCORE]]
            cols1 = pd.MultiIndex.from_product(iterable_1,
                                               names=['Samples', 'Stats'])
            iterable_2 = [['Details'], ['Rank', 'Name']]
            cols2 = pd.MultiIndex.from_product(iterable_2)
            cols = cols1.append(cols2)
            data_frame = pd.DataFrame.from_items(list_rows,
                                                 orient='index',
                                                 columns=cols)
            data_frame.index.names = ['TaxId']
            data_frame.to_excel(xlsxwriter, sheet_name=str(excel))
        elif excel is Excel.CMPLXCRUNCHER:
            target_ranks: List = [Rank.NO_RANK]
            if args.controls:
                target_ranks = [Rank.SPECIES, Rank.GENUS,  # Ranks of interest
                                Rank.FAMILY, Rank.ORDER]  # for cmplxcruncher
            for rank in target_ranks:  # Once for no rank dependency (NO_RANK)
                indexes: List[int]
                sheet_name: str
                columns: List[str]
                if args.controls:
                    indexes = [i for i in range(len(raw_samples), len(samples))
                               if (samples[i].startswith(STR_CONTROL)
                                   and rank.name.lower() in samples[i])]
                    sheet_name = f'{STR_CONTROL}_{rank.name.lower()}'
                    columns = [samples[i].split('_')[2] for i in indexes]
                else:  # No rank dependency
                    indexes = list(range(len(raw_samples)))
                    sheet_name = f'raw_samples_{rank.name.lower()}'
                    columns = [samples[i].split('_')[0] for i in indexes]
                list_rows = []
                polytree.to_items(taxonomy=ncbi, items=list_rows,
                                  sample_indexes=indexes)
                data_frame = pd.DataFrame.from_items(list_rows,
                                                     orient='index',
                                                     columns=columns)
                data_frame.index.names = ['TaxId']
                data_frame.to_excel(xlsxwriter, sheet_name=sheet_name)
        else:
            raise Exception(red('\nERROR!'),
                            f'Unknown Excel option "{excel}"')
        xlsxwriter.save()
        print(green('OK!'))

    # timing initialization
    start_time: float = time.time()
    # Program header
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Parse arguments
    argparser = configure_parser()
    args = argparser.parse_args()
    outputs: List[Filename] = args.file
    reports: List[Filename] = args.report
    lmats: List[Filename] = args.lmat
    input_files: List[Filename]
    nodesfile: Filename = Filename(os.path.join(args.nodespath, NODES_FILE))
    namesfile: Filename = Filename(os.path.join(args.nodespath, NAMES_FILE))
    htmlfile: Filename = args.outhtml
    collapse: bool = not args.nokollapse
    excluding: Set[TaxId] = set(args.exclude)
    including: Set[TaxId] = set(args.include)
    scoring: Scoring = Scoring[args.scoring]
    excel: Excel = Excel[args.excel]

    check_debug()

    plasmidfile: Filename = None
    process: Callable[..., Tuple[Sample, TaxTree, SampleDataByTaxId, Err]]
    select_inputs()
    check_controls()
    if not htmlfile:
        select_html_file()

    # Load NCBI nodes, names and build children
    ncbi: Taxonomy = Taxonomy(nodesfile, namesfile, plasmidfile,
                              collapse, excluding, including, args.debug)

    # If dummy flag enabled, just create dummy krona and exit
    if args.dummy:
        _debug_dummy_plot(ncbi, htmlfile, scoring)
        exit(0)

    # Declare variables that will hold results for the samples analyzed
    trees: Dict[Sample, TaxTree] = {}
    counts: Dict[Sample, Counter[TaxId]] = {}
    accs: Dict[Sample, Counter[TaxId]] = {}
    taxids: Dict[Sample, TaxLevels] = {}
    scores: Dict[Sample, Dict[TaxId, Score]] = {}
    samples: List[Sample] = []
    raw_samples: List[Sample] = []

    # Define dictionary of parameters for methods to be called (to be extended)
    kwargs = {'controls': args.controls,
              'ctrlminscore': (
                  args.ctrlminscore
                  if args.ctrlminscore is not None else args.minscore),
              'ctrlmintaxa': (
                  args.ctrlmintaxa
                  if args.ctrlmintaxa is not None else args.mintaxa),
              'debug': args.debug, 'root':args.takeoutroot,
              'lmat': bool(lmats), 'minscore': args.minscore,
              'mintaxa': args.mintaxa, 'scoring': scoring, 'taxonomy': ncbi,
              }
    # The big stuff (done in parallel)
    read_samples()
    # Avoid cross analysis if just one report file or explicitly stated by flag
    if len(raw_samples) > 1 and not args.avoidcross:
        analyze_samples()
        if args.summary != 'avoid':
            summaries: List[Sample] = []
            summarize_samples()
            if args.summary == 'only':
                samples = raw_samples + summaries
            else:
                samples.extend(summaries)
    # Final result generation is done in sequential mode

    polytree: MultiTree = MultiTree(samples=samples)
    generate_krona()
    if _USE_PANDAS:
        generate_excel()
    else:
        print(yellow('WARNING!'),
              'Pandas not installed: Excel cannot be created.')

    # Timing results
    print(gray('Total elapsed time:'), time.strftime(
        "%H:%M:%S", time.gmtime(time.time()-start_time)))


if __name__ == '__main__':
    main()
