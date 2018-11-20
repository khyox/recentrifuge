#!/usr/bin/env python3
#
#     Copyright (C) 2017, 2018, Jose Manuel Martí Martínez
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program. If not, see <https://www.gnu.org/licenses/>.
#
"""
Analyze metagenomic taxonomic classification data.
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

from recentrifuge.centrifuge import select_centrifuge_inputs
from recentrifuge.clark import select_clark_inputs
from recentrifuge.config import Filename, Sample, Id, Score, Scoring, Extra
from recentrifuge.config import HTML_SUFFIX, DEFMINTAXA, TAXDUMP_PATH
from recentrifuge.config import LICENSE, Err, Classifier
from recentrifuge.config import NODES_FILE, NAMES_FILE, PLASMID_FILE
from recentrifuge.config import STR_CONTROL, STR_EXCLUSIVE, STR_SHARED
from recentrifuge.config import STR_CONTROL_SHARED
from recentrifuge.config import gray, red, green, yellow, blue, magenta
from recentrifuge.core import process_rank, summarize_analysis
from recentrifuge.krona import COUNT, UNASSIGNED, SCORE
from recentrifuge.krona import KronaTree
from recentrifuge.lmat import select_lmat_inputs
from recentrifuge.rank import Rank, TaxLevels
from recentrifuge.stats import SampleStats
from recentrifuge.taxclass import process_output, process_report
from recentrifuge.taxonomy import Taxonomy
from recentrifuge.trees import TaxTree, MultiTree, SampleDataById

# optional package pandas (to generate extra output)
_USE_PANDAS = True
try:
    import pandas as pd
except ImportError:
    pd = None
    _USE_PANDAS = False

__version__ = '0.22.3'
__author__ = 'Jose Manuel Martí'
__date__ = 'Nov 2018'


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
    polytree.grow(ontology=taxonomy)
    polytree.toxml(ontology=taxonomy, krona=krona)
    krona.tohtml(htmlfile, pretty=True)
    print(green('OK!'))


def main():
    """Main entry point to Recentrifuge."""

    def vprint(*arguments) -> None:
        """Print only if verbose/debug mode is enabled"""
        if args.debug:
            print(*arguments)

    def configure_parser():
        """Argument Parser Configuration"""
        parser = argparse.ArgumentParser(
            description='Analyze results of metagenomic taxonomic classifiers',
            epilog=f'%(prog)s  - Release {__version__} - {__date__}' + LICENSE,
            formatter_class=argparse.RawDescriptionHelpFormatter
        )
        parser.add_argument(
            '-V', '--version',
            action='version',
            version=f'%(prog)s version {__version__} released in {__date__}'
        )
        parser_in = parser.add_argument_group(
            'input', 'Define Recentrifuge input files and formats')
        parser_in.add_argument(
            '-n', '--nodespath',
            action='store',
            metavar='PATH',
            default=TAXDUMP_PATH,
            help=('path for the nodes information files '
                  '(nodes.dmp and names.dmp from NCBI)')
        )
        parser_filein = parser_in.add_mutually_exclusive_group(required=True)
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
        parser_filein.add_argument(
            '-k', '--clark',
            action='append',
            metavar='FILE',
            type=Filename,
            help=('CLARK(S) output files. If a single directory is entered, '
                  'every .csv file inside will be taken as a different sample.'
                  ' Multiple -k is available to include several samples.')
        )
        parser_filein.add_argument(
            '-r', '--report',
            action='append',
            metavar='FILE',
            type=Filename,
            help=('Centrifuge/Kraken report files '
                  '(multiple -r is available to include several samples)')
        )
        parser_out = parser.add_argument_group(
            'output', 'Related to the Recentrifuge output files')
        parser_out.add_argument(
            '-o', '--outhtml',
            action='store',
            metavar='FILE',
            type=Filename,
            help='HTML output file (if not given, the filename will be '
                 'inferred from input files)'
        )
        parser_out.add_argument(
            '-e', '--extra',
            action='store',
            metavar='OUTPUT_TYPE',
            choices=[str(extra) for extra in Extra],
            default=str(Extra(0)),
            help=(f'type of extra output to be generated, and can be one of '
                  f'{[str(extra) for extra in Extra]}')
        )
        parser_coarse = parser.add_argument_group(
            'tuning',
            'Coarse tuning of algorithm parameters')
        parser_cross = parser_coarse.add_mutually_exclusive_group(
            required=False)
        parser_cross.add_argument(
            '-c', '--controls',
            action='store',
            metavar='CONTROLS_NUMBER',
            type=int,
            default=0,
            help=('this number of first samples will be treated as negative '
                  'controls; default is no controls')
        )
        parser_coarse.add_argument(
            '-s', '--scoring',
            action='store',
            metavar='SCORING',
            choices=[str(each_score) for each_score in Scoring],
            default=str(Scoring(0)),
            help=(f'type of scoring to be applied, and can be one of '
                  f'{[str(scoring) for scoring in Scoring]}')
        )
        parser_coarse.add_argument(
            '-y', '--minscore',
            action='store',
            metavar='NUMBER',
            type=lambda txt: Score(float(txt)),
            default=None,
            help=('minimum score/confidence of the classification of a read '
                  'to pass the quality filter; all pass by default')
        )
        parser_coarse.add_argument(
            '-m', '--mintaxa',
            action='store',
            metavar='INT',
            type=int,
            default=DEFMINTAXA,
            help='minimum taxa to avoid collapsing one level to the parent one'
        )
        parser_coarse.add_argument(
            '-x', '--exclude',
            action='append',
            metavar='TAXID',
            type=Id,
            default=[],
            help=('NCBI taxid code to exclude a taxon and all underneath '
                  '(multiple -x is available to exclude several taxid)')
        )
        parser_coarse.add_argument(
            '-i', '--include',
            action='append',
            metavar='TAXID',
            type=Id,
            default=[],
            help=('NCBI taxid code to include a taxon and all underneath '
                  '(multiple -i is available to include several taxid); '
                  'by default, all the taxa are considered for inclusion')
        )
        parser_cross.add_argument(
            '-a', '--avoidcross',
            action='store_true',
            help='avoid cross analysis'
        )
        parser_fine = parser.add_argument_group(
            'fine tuning',
            'Fine tuning of algorithm parameters')
        parser_fine.add_argument(
            '-z', '--ctrlminscore',
            action='store',
            metavar='NUMBER',
            type=lambda txt: Score(float(txt)),
            default=None,
            help=('minimum score/confidence of the classification of a read '
                  'in control samples to pass the quality filter; if defaults '
                  'to "minscore"')
        )
        parser_fine.add_argument(
            '-w', '--ctrlmintaxa',
            action='store',
            metavar='INT',
            type=int,
            default=None,
            help='minimum taxa to avoid collapsing one level to the parent one'
                 ' in control samples; it defaults to "mintaxa"'
        )
        parser_fine.add_argument(
            '-u', '--summary',
            action='store',
            metavar='OPTION',
            choices=['add', 'only', 'avoid'],
            default='add',
            help=('select to "add" summary samples to other samples, or to '
                  '"only" show summary samples or to "avoid" summaries at all')
        )
        parser_fine.add_argument(
            '-t', '--takeoutroot',
            action='store_true',
            help='remove counts directly assigned to the "root" level'
        )
        parser_fine.add_argument(
            '--nokollapse',
            action='store_true',
            help='show the "cellular organisms" taxon'
        )
        parser_mode = parser.add_argument_group('advanced',
                                                'Advanced modes of running')
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
        return parser

    def check_debug():
        """Check debugging mode"""
        if args.debug:
            print(blue('INFO:'), gray('Debugging mode activated'))
            print(blue('INFO:'), gray('Active parameters:'))
            for key, value in vars(args).items():
                if value:
                    print(gray(f'\t{key} ='), f'{value}')

    def select_inputs():
        """Choose right classifier, input and output files"""
        nonlocal process, scoring, input_files, plasmidfile, classifier

        if reports:
            classifier = Classifier.KRAKEN
            process = process_report
            input_files = reports
        elif clarks:
            classifier = Classifier.CLARK
            process = process_output
            input_files = clarks
            if len(clarks) == 1 and os.path.isdir(clarks[0]):
                select_clark_inputs(clarks)
        elif lmats:
            classifier = Classifier.LMAT
            scoring = Scoring.LMAT
            process = process_output
            input_files = lmats
            plasmidfile = Filename(os.path.join(args.nodespath, PLASMID_FILE))
            select_lmat_inputs(lmats)
        elif outputs:
            classifier = Classifier.CENTRIFUGE
            process = process_output
            input_files = outputs
            if len(outputs) == 1 and os.path.isdir(outputs[0]):
                select_centrifuge_inputs(outputs)

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
                for file, (sample, tree, out, stat, err) in zip(
                        input_files, [r.get() for r in async_results]):
                    if err is Err.NO_ERROR:
                        samples.append(sample)
                        trees[sample] = tree
                        taxids[sample] = out.get_taxlevels()
                        counts[sample] = out.counts
                        accs[sample] = out.accs
                        scores[sample] = out.scores
                        stats[sample] = stat
                    elif err is Err.VOID_CTRL:
                        print('There were void controls.', red('Aborting!'))
                        exit(1)
        else:  # sequential processing of each sample
            for num, file in enumerate(input_files):
                (sample, tree, out, stat, err) = process(
                    file, True if num < args.controls else False, **kwargs)
                if err is Err.NO_ERROR:
                    samples.append(sample)
                    trees[sample] = tree
                    taxids[sample] = out.get_taxlevels()
                    counts[sample] = out.counts
                    accs[sample] = out.accs
                    scores[sample] = out.scores
                    stats[sample] = stat
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
                                     )
        polytree.grow(ontology=ncbi,
                      abundances=counts,
                      accs=accs,
                      scores=scores)
        print(green('OK!'))
        print(gray('Generating final plot (') + magenta(htmlfile) +
              gray(')... '), end='')
        sys.stdout.flush()
        polytree.toxml(ontology=ncbi, krona=krona)
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

        # Save basic statistics of raw samples
        data_frame: pd.DataFrame = pd.DataFrame.from_dict(
            {raw: stats[raw].to_dict() for raw in raw_samples})
        if extra is Extra.FULL or extra is Extra.CMPLXCRUNCHER:
            data_frame.to_excel(xlsxwriter, sheet_name='_sample_stats')
        elif extra is Extra.CSV or extra is Extra.TSV:
            data_frame.to_csv(sv_base + '.stat.' + sv_ext, **sv_kargs)

        # Save taxid related statistics per sample
        if extra is Extra.FULL or extra is Extra.CSV or extra is Extra.TSV:
            polytree.to_items(ontology=ncbi, items=list_rows)
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
                target_ranks.extend(Rank.selected_ranks)
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
                polytree.to_items(ontology=ncbi, items=list_rows,
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
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} - {__date__}'
          f' =-= by {__author__} =-=\n')
    sys.stdout.flush()

    # Parse arguments
    argparser = configure_parser()
    args = argparser.parse_args()
    outputs: List[Filename] = args.file
    reports: List[Filename] = args.report
    lmats: List[Filename] = args.lmat
    clarks: List[Filename] = args.clark
    input_files: List[Filename]
    nodesfile: Filename = Filename(os.path.join(args.nodespath, NODES_FILE))
    namesfile: Filename = Filename(os.path.join(args.nodespath, NAMES_FILE))
    htmlfile: Filename = args.outhtml
    collapse: bool = not args.nokollapse
    excluding: Set[Id] = set(args.exclude)
    including: Set[Id] = set(args.include)
    scoring: Scoring = Scoring[args.scoring]
    extra: Extra = Extra[args.extra]

    check_debug()

    plasmidfile: Filename = None
    classifier: Classifier
    process: Callable[..., Tuple[Sample, TaxTree, SampleDataById,
                                 SampleStats, Err]]
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
    counts: Dict[Sample, Counter[Id]] = {}
    accs: Dict[Sample, Counter[Id]] = {}
    taxids: Dict[Sample, TaxLevels] = {}
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
                  if args.ctrlmintaxa is not None else args.mintaxa),
              'debug': args.debug, 'root': args.takeoutroot,
              'classifier': classifier, 'minscore': args.minscore,
              'mintaxa': args.mintaxa, 'scoring': scoring, 'ontology': ncbi,
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
              'Pandas not installed: Extra output cannot be saved.')

    # Timing results
    print(gray('Total elapsed time:'), time.strftime(
        "%H:%M:%S", time.gmtime(time.time() - start_time)))


if __name__ == '__main__':
    main()
