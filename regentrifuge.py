#!/usr/bin/env python3
"""
Post-process Gene Ontology output.
"""
# pylint: disable=no-name-in-module, not-an-iterable
import argparse
import collections as col
import multiprocessing as mp
import os
import platform
import sys
import time
from typing import Counter, List, Dict, Set

from recentrifuge.config import DEFMINGENE, GODUMP_PATH, NODES_FILE, NAMES_FILE
from recentrifuge.config import Filename, Sample, Id, Score, Scoring, Chart
from recentrifuge.config import HTML_SUFFIX, Extra, Err
from recentrifuge.config import STR_CONTROL, STR_EXCLUSIVE, STR_SHARED
from recentrifuge.config import STR_CONTROL_SHARED
from recentrifuge.config import gray, red, green, yellow, blue, magenta
from recentrifuge.core import process_rank, summarize_analysis
from recentrifuge.gene_ontology import GeneOntology
from recentrifuge.go import process_goutput
from recentrifuge.krona import KronaTree
from recentrifuge.krona import COUNT, UNASSIGNED, SCORE
from recentrifuge.rank import Rank, TaxLevels
from recentrifuge.stats import SampleStats
from recentrifuge.trees import TaxTree, MultiTree

# optional package pandas (to generate Extra output)
_USE_PANDAS = True
try:
    import pandas as pd
except ImportError:
    pd = None
    _USE_PANDAS = False

__version__ = '0.2.0'
__author__ = 'Jose Manuel Marti'
__date__ = 'Nov 2018'


def _debug_dummy_plot(geno: GeneOntology,
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
    polytree.grow(ontology=geno)
    polytree.toxml(ontology=geno, krona=krona)
    krona.tohtml(htmlfile, pretty=True)
    print(green('OK!'))


def main():
    """Main entry point to Regentrifuge."""

    def vprint(*arguments) -> None:
        """Print only if verbose/debug mode is enabled"""
        if args.debug:
            print(*arguments)

    def configure_parser():
        """Argument Parser Configuration"""
        parser = argparse.ArgumentParser(
            description='Post-process Gene Ontology output',
            epilog=f'%(prog)s  - {__author__} - {__date__}',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
        parser.add_argument(
            '-n', '--nodespath',
            action='store',
            metavar='PATH',
            default=GODUMP_PATH,
            help=('path for the nodes information files '
                  '(nodes.dmp and names.dmp adapted from GO)')
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
            help=('GO output files in tab separated columns with the format:'
                  'taxid, gi, go, score and evalue. Currently, the data for '
                  'the analysis is retrieved from the last 3 columns. '
                  ' Multiple -f is available to include several samples.')
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
                  'controls; default is no controls. CAUTION! Regentrifuge '
                  'direct support of control samples is experimental.')
        )
        parser_output = parser.add_argument_group(
            'output', 'Related to the output files')
        parser_output.add_argument(
            '-e', '--extra',
            action='store',
            metavar='OUTPUT_TYPE',
            choices=[str(extra) for extra in Extra],
            default=str(Extra(0)),
            help=(f'type of scoring to be applied, and can be one of '
                  f'{[str(extra) for extra in Extra]}')
        )
        parser_output.add_argument(
            '-o', '--outhtml',
            action='store',
            metavar='FILE',
            type=Filename,
            help='HTML output file (if not given the filename will be '
                 'inferred from input files)',
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
            '-m', '--mingene',
            action='store',
            metavar='INT',
            type=int,
            default=DEFMINGENE,
            help='minimum genes to avoid collapsing GO genes hierarchy levels'
        )
        parser_tuning.add_argument(
            '-i', '--include',
            action='append',
            metavar='TAXID',
            type=Id,
            default=[],
            help=('GO code to include a GO node and all underneath '
                  '(multiple -i is available to include several GOs); '
                  'by default all the GOs are considered for inclusion')
        )
        parser_tuning.add_argument(
            '-x', '--exclude',
            action='append',
            metavar='TAXID',
            type=Id,
            default=[],
            help=('GO code to exclude a GO node and all underneath '
                  '(multiple -x is available to exclude several GOs)')
        )
        parser_tuning.add_argument(
            '-s', '--scoring',
            action='store',
            metavar='SCORING',
            choices=[str(each_score) for each_score in Scoring],
            default=str(Scoring(0)),
            help=argparse.SUPPRESS,  # Different scoring still not supported
            # help=(f'type of scoring to be applied, and can be one of '
            #      f'{[str(scoring) for scoring in Scoring]}')
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
            '-t', '--takeoutroot',
            action='store_true',
            help='remove counts directly assigned to the "root" level'
        )
        parser_tuning.add_argument(
            '-y', '--minscore',
            action='store',
            metavar='NUMBER',
            type=lambda txt: Score(float(txt)),
            default=None,
            help=('minimum score/confidence of the annotation of a scaffold '
                  'to pass the quality filter; all pass by default')
        )
        return parser

    def check_debug():
        """Check debugging mode"""
        if args.debug:
            print(gray('INFO: Debugging mode activated\n'))

    def read_samples():
        """Read samples"""
        print(gray('\nPlease, wait, processing files in parallel...\n'))
        # Enable parallelization with 'spawn' under known platforms
        if platform.system() and not args.sequential:  # Only for known systems
            mpctx = mp.get_context('fork')
            with mpctx.Pool(processes=min(os.cpu_count(),
                                          len(input_files))) as pool:
                async_results = [pool.apply_async(
                    process_goutput,
                    args=[input_files[num],  # file name
                          True if num < args.controls else False],  # is ctrl?
                    kwds=kwargs
                ) for num in range(len(input_files))]
                for file, (sample, tree, out, stat, err) in zip(
                        input_files, [r.get() for r in async_results]):
                    if err is Err.NO_ERROR:
                        samples.append(sample)
                        trees[sample] = tree
                        ids[sample] = out.get_taxlevels()
                        counts[sample] = out.counts
                        accs[sample] = out.accs
                        scores[sample] = out.scores
                        stats[sample] = stat
                    elif err is Err.VOID_CTRL:
                        print('There were void controls.', red('Aborting!'))
                        exit(1)
        else:  # sequential processing of each sample
            for num, file in enumerate(input_files):
                (sample, tree, out, stat, err) = process_goutput(
                    file, True if num < args.controls else False, **kwargs)
                if err is Err.NO_ERROR:
                    samples.append(sample)
                    trees[sample] = tree
                    ids[sample] = out.get_taxlevels()
                    counts[sample] = out.counts
                    accs[sample] = out.accs
                    scores[sample] = out.scores
                    stats[sample] = stat
                elif err is Err.VOID_CTRL:
                    print('There were void controls.', red('Aborting!'))
                    exit(1)
        raw_samples.extend(samples)  # Store raw sample names

    def analyze_samples():
        """Cross analysis of samples in parallel by rank"""
        print(gray('Please, wait. Performing cross analysis in parallel...\n'))
        # Update kwargs with more parameters for the followings func calls
        kwargs.update({'taxids': ids, 'counts': counts, 'scores': scores,
                       'accs': accs, 'raw_samples': raw_samples})
        if platform.system() and not args.sequential:  # Only for known systems
            mpctx = mp.get_context('fork')  # Important for OSX&Win
            with mpctx.Pool(processes=min(os.cpu_count(), len(
                    Rank.genomic_ranks))) as pool:
                async_results = [pool.apply_async(
                    process_rank,
                    args=[level],
                    kwds=kwargs
                ) for level in Rank.genomic_ranks]
                for level, (smpls, abunds, accumulators, score) in zip(
                        Rank.genomic_ranks,
                        [r.get() for r in async_results]):
                    samples.extend(smpls)
                    counts.update(abunds)
                    accs.update(accumulators)
                    scores.update(score)
        else:  # sequential processing of each selected rank
            for level in Rank.genomic_ranks:
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

        print(gray('\nBuilding the ontology multiple tree... '), end='')
        sys.stdout.flush()
        krona: KronaTree = KronaTree(samples,
                                     num_raw_samples=len(raw_samples),
                                     stats=stats,
                                     min_score=Score(
                                         min([min(scores[sample].values())
                                              for sample in samples
                                              if len(scores[sample])])),
                                     max_score=Score(
                                         max([max(scores[sample].values())
                                              for sample in samples
                                              if len(scores[sample])])),
                                     scoring=scoring,
                                     chart=Chart.GENOMIC,
                                     )
        polytree.grow(ontology=geno,
                      abundances=counts,
                      accs=accs,
                      scores=scores)
        print(green('OK!'))
        print(gray('Generating final plot (') + magenta(htmlfile) +
              gray(')... '), end='')
        sys.stdout.flush()
        polytree.toxml(ontology=geno, krona=krona)
        krona.tohtml(htmlfile, pretty=False)
        print(green('OK!'))

    def save_extra_output():
        """Save extra output with results via pandas DataFrame"""

        # Initial check and info
        if extra is Extra.FULL or extra is Extra.CMPLXCRUNCHER:
            vprint(blue('INFO:'),
                   gray('Saving extra output as an Excel file.'))
        elif extra is Extra.CSV:
            vprint(blue('INFO:'),
                   gray('Saving extra output as CSV files.'))
        elif extra is Extra.TSV:
            vprint(blue('INFO:'),
                   gray('Saving extra output as TSV files.'))
        else:
            raise Exception(red('\nERROR!'), f'Unknown Extra option "{extra}"')

        # Setup of extra file name
        xlsxwriter: pd.ExcelWriter = None
        if extra is Extra.FULL or extra is Extra.CMPLXCRUNCHER:
            xlsx_name: Filename = Filename(htmlfile.split('.html')[0] +
                                           '.xlsx')
            print(gray(f'Generating Excel {str(extra).lower()} summary (') +
                  magenta(xlsx_name) + gray(')... '), end='')
            sys.stdout.flush()
            xlsxwriter = pd.ExcelWriter(xlsx_name)
        elif extra is Extra.CSV or Extra.TSV:
            sv_base: Filename = Filename(htmlfile.split('.html')[0])
            sv_ext: Filename
            sv_kargs: Dict[str, str]
            if extra is Extra.CSV:
                sv_ext = 'csv'
                sv_kargs = {}
            elif extra is Extra.TSV:
                sv_ext = 'tsv'
                sv_kargs = {'sep': '\t'}
            print(gray(f'Generating {str(extra).lower()} extra output (') +
                  magenta(sv_base + '.*.' + sv_ext) + gray(')... '), end='')
            sys.stdout.flush()
        list_rows: List = []

        # Save raw samples basic statistics
        data_frame: pd.DataFrame = pd.DataFrame.from_dict(
            {raw: stats[raw].to_dict() for raw in raw_samples})
        if extra is Extra.FULL or extra is Extra.CMPLXCRUNCHER:
            data_frame.to_excel(xlsxwriter, sheet_name='_sample_stats')
        elif extra is Extra.CSV or extra is Extra.TSV:
            data_frame.to_csv(sv_base + '.stat.' + sv_ext, **sv_kargs)

        # Save taxid related statistics per sample
        if extra is Extra.FULL or extra is Extra.CSV or extra is Extra.TSV:
            polytree.to_items(ontology=geno, items=list_rows)
            # Generate the pandas DataFrame from items and export to Extra
            iterable_1 = [samples, [COUNT, UNASSIGNED, SCORE]]
            cols1 = pd.MultiIndex.from_product(iterable_1,
                                               names=['Samples', 'Stats'])
            iterable_2 = [['Details'], ['Rank', 'Name']]
            cols2 = pd.MultiIndex.from_product(iterable_2)
            cols = cols1.append(cols2)
            data_frame = pd.DataFrame.from_items(list_rows,
                                                 orient='index',
                                                 columns=cols)
            data_frame.index.names = ['Id']
            if extra is Extra.CSV or extra is Extra.TSV:
                data_frame.to_csv(sv_base + '.data.' + sv_ext, **sv_kargs)
            else:
                data_frame.to_excel(xlsxwriter, sheet_name=str(extra))

        elif extra is Extra.CMPLXCRUNCHER:
            target_ranks: List = [Rank.NO_RANK]
            if args.controls:  # if controls, add specific sheet for rank
                target_ranks.extend(Rank.genomic_ranks)
            for rank in target_ranks:  # Once for no rank dependency (NO_RANK)
                indexes: List[int]
                sheet_name: str
                columns: List[str] = ['__Rank', '__Name']
                if args.controls:
                    indexes = [i for i in range(len(raw_samples), len(samples))
                               # Check if sample ends in _(STR_CONTROL)_(rank)
                               if (STR_CONTROL in samples[i].split('_')[-2:]
                                   and rank.name.lower() in
                                   samples[i].split('_')[-1:])]
                    sheet_name = f'{STR_CONTROL}_{rank.name.lower()}'
                    columns.extend([samples[i].replace(
                        '_' + STR_CONTROL + '_' + rank.name.lower(), '')
                        for i in indexes])
                if rank is Rank.NO_RANK:  # No rank dependency
                    indexes = list(range(len(raw_samples)))
                    sheet_name = f'raw_samples_{rank.name.lower()}'
                    columns.extend(raw_samples)
                list_rows = []
                polytree.to_items(ontology=geno, items=list_rows,
                                  sample_indexes=indexes)
                data_frame = pd.DataFrame.from_items(list_rows,
                                                     orient='index',
                                                     columns=columns)
                data_frame.index.names = ['Id']
                data_frame.to_excel(xlsxwriter, sheet_name=sheet_name)

        else:
            raise Exception(red('\nERROR!'), f'Unknown Extra option "{extra}"')
        if xlsxwriter is not None:
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
    input_files: List[Filename]
    nodesfile: Filename = Filename(os.path.join(args.nodespath, NODES_FILE))
    namesfile: Filename = Filename(os.path.join(args.nodespath, NAMES_FILE))
    htmlfile: Filename = args.outhtml
    scoring: Scoring = Scoring[args.scoring]
    args.ctrlminscore = None
    args.ctrlmingene = None
    excluding: Set[Id] = set(args.exclude)
    including: Set[Id] = set(args.include)
    extra: Extra = Extra[args.extra]

    check_debug()
    input_files = outputs  # Done by select_inputs() in recentrifuge
    if htmlfile is None:   # Done by select_html_file() in recentrifuge
        htmlfile = Filename(outputs[0].split('_go')[0] + HTML_SUFFIX)

    # Load GO nodes, names and build children
    geno: GeneOntology = GeneOntology(nodesfile, namesfile,
                                    excluding, including, args.debug)

    # If dummy flag enabled, just create dummy krona and exit
    if args.dummy:
        _debug_dummy_plot(geno, htmlfile, scoring)
        exit(0)

    # Declare variables that will hold results for the samples analyzed
    trees: Dict[Sample, TaxTree] = {}
    counts: Dict[Sample, Counter[Id]] = {}
    accs: Dict[Sample, Counter[Id]] = {}
    ids: Dict[Sample, TaxLevels] = {}
    scores: Dict[Sample, Dict[Id, Score]] = {}
    stats: Dict[Sample, SampleStats] = {}
    samples: List[Sample] = []
    raw_samples: List[Sample] = []

    # Define dictionary of parameters for methods to be called (to be extended)
    kwargs = {'controls': args.controls,
              'ctrlminscore': (
                  args.ctrlminscore
                  if args.ctrlminscore is not None else args.minscore),
              'ctrlmintaxa': (
                  args.ctrlmintaxa
                  if args.ctrlmingene is not None else args.mingene),
              'debug': args.debug, 'root': args.takeoutroot,
              'minscore': args.minscore,
              'mintaxa': args.mingene, 'scoring': scoring, 'ontology': geno,
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
        save_extra_output()
    else:
        print(yellow('WARNING!'),
              'Pandas not installed: Extra cannot be created.')

    # Timing results
    print(gray('Total elapsed time:'), time.strftime(
        "%H:%M:%S", time.gmtime(time.time() - start_time)))


if __name__ == '__main__':
    main()
