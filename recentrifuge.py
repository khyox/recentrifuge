#!/usr/bin/env python3
"""
Post-process Centrifuge/Kraken output.
"""
# pylint: disable=no-name-in-module, not-an-iterable
import argparse
import multiprocessing as mp
import os
import platform
import sys
from typing import Counter, List, Dict, Set, Callable, Optional, Tuple

# optional package pandas (to generate Excel output)
_use_pandas = True
try:
    import pandas as pd
except ImportError:
    _use_pandas = False

from recentrifuge.centrifuge import process_report, process_output
from recentrifuge.config import Filename, Sample, TaxId, Score, Scoring, Excel
from recentrifuge.config import NODES_FILE, NAMES_FILE, TAXDUMP_PATH
from recentrifuge.config import HTML_SUFFIX, DEFMINTAXA, STR_CONTROL
from recentrifuge.core import Taxonomy, TaxLevels, TaxTree, MultiTree, Rank
from recentrifuge.core import process_rank
from recentrifuge.krona import KronaTree
from recentrifuge.krona import COUNT, UNASSIGNED, SCORE

__version__ = '0.12.4'
__author__ = 'Jose Manuel Marti'
__date__ = 'Nov 2017'


def _debug_dummy_plot(taxonomy: Taxonomy,
                      htmlfile: Filename,
                      scoring: Scoring = Scoring.SHEL,
                      ):
    """

    Generate dummy Krona plot via Krona 2.0 XML spec and finish

    """
    print(f'\033[90mGenerating dummy Krona plot {htmlfile}...\033[0m', end='')
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
    print('\033[92m OK! \033[0m')
    # End execution
    quit()


def main():
    """Main entry point to recentrifuge."""
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Post-process Centrifuge/Kraken output',
        epilog=f'%(prog)s  - {__author__} - {__date__}',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s release {__version__} ({__date__})'
    )
    filein = parser.add_mutually_exclusive_group(required=True)
    filein.add_argument(
        '-f', '--file',
        action='append',
        metavar='FILE',
        help=('Centrifuge output files '
              '(multiple -f is available to include several samples in plot)')
    )
    filein.add_argument(
        '-r', '--report',
        action='append',
        metavar='FILE',
        help=('Centrifuge/Kraken report files '
              '(multiple -r is available to include several samples in plot)')
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        help='increase output verbosity'
    )
    parser.add_argument(
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default=TAXDUMP_PATH,
        help=('path for the nodes information files (nodes.dmp and names.dmp'
              ' from NCBI)')
    )
    parser.add_argument(
        '-m', '--mintaxa',
        action='store',
        metavar='INT',
        default=DEFMINTAXA,
        help=('minimum taxa to avoid collapsing one level to the parent one')
    )
    parser.add_argument(
        '-k', '--nokollapse',
        action='store_true',
        help='show the "cellular organisms" taxon'
    )
    parser.add_argument(
        '-o', '--outhtml',
        action='store',
        metavar='FILE',
        default=None,
        help='HTML output file (if not given the filename will be '
             'inferred from input files)'
    )
    parser.add_argument(
        '-i', '--include',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to include a taxon and all underneath '
              '(multiple -i is available to include several taxid). '
              'By default all the taxa is considered for inclusion.')
    )
    parser.add_argument(
        '-x', '--exclude',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to exclude a taxon and all underneath '
              '(multiple -x is available to exclude several taxid)')
    )
    parser.add_argument(
        '-s', '--scoring',
        action='store',
        metavar='SCORING',
        choices=[str(scoring) for scoring in Scoring],
        default=str(Scoring(0)),
        help=(f'type of scoring to be applied, and can be one of '
              f'{[str(scoring) for scoring in Scoring]}')
    )
    parser.add_argument(
        '--sequential',
        action='store_true',
        help='deactivate parallel processing'
    )
    parser.add_argument(
        '-a', '--avoidcross',
        action='store_true',
        help='avoid cross analysis'
    )
    parser.add_argument(
        '-c', '--control',
        action='store_true',
        help='take the first sample as negative control'
    )
    parser.add_argument(
        '-e', '--excel',
        action='store',
        metavar='OUTPUT_TYPE',
        choices=[str(excel) for excel in Excel],
        default=str(Excel(0)),
        help=(f'type of scoring to be applied, and can be one of '
              f'{[str(excel) for excel in Excel]}')
    )

    # Parse arguments
    args = parser.parse_args()
    outputs = args.file
    files = args.report
    verb = args.verbose
    nodesfile = os.path.join(args.nodespath, NODES_FILE)
    namesfile = os.path.join(args.nodespath, NAMES_FILE)
    mintaxa = int(args.mintaxa)
    collapse = not args.nokollapse
    excluding: Set[TaxId] = set(args.exclude)
    including: Set[TaxId] = set(args.include)
    scoring: Scoring = Scoring[args.scoring]
    sequential = args.sequential
    avoidcross = args.avoidcross
    control = args.control
    excel: Excel = Excel[args.excel]
    htmlfile: Filename = args.outhtml
    if not htmlfile:
        if files:
            htmlfile = files[0].split('_mhl')[0] + HTML_SUFFIX
        else:
            htmlfile = outputs[0].split('_mhl')[0] + HTML_SUFFIX

    # Program header
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Load NCBI nodes, names and build children
    ncbi: Taxonomy = Taxonomy(nodesfile, namesfile,
                              collapse, excluding, including)

    #  _debug_dummy_plot(ncbi, htmlfile, scoring)

    # Declare variables that will hold results for the samples analyzed
    trees: Dict[Sample, TaxTree] = {}
    abundances: Dict[Sample, Counter[TaxId]] = {}
    accs: Dict[Sample, Counter[TaxId]] = {}
    taxids: Dict[Sample, TaxLevels] = {}
    scores: Dict[Sample, Dict[TaxId, Score]] = {}
    samples: List[Sample] = []
    raw_samples: List[Sample] = []
    #
    # Processing of input files in parallel
    #
    process: Callable[..., Tuple[Sample, TaxTree, TaxLevels,
                                 Counter[TaxId], Counter[TaxId],
                                 Dict[TaxId, Optional[Score]]]]
    if files:
        process = process_report
        files = files
    else:
        process = process_output
        files = outputs

    print('\033[90mPlease, wait, processing files in parallel...\033[0m\n')
    kwargs = {'taxonomy': ncbi, 'mintaxa': mintaxa,
              'verb': verb, 'scoring': scoring}
    # Enable parallelization with 'spawn' under known platforms
    if platform.system() and not sequential:  # Only for known platforms
        mpctx = mp.get_context('spawn')  # Important for OSX&Win
        with mpctx.Pool(processes=min(os.cpu_count(),
                                      len(files))) as pool:
            async_results = [pool.apply_async(
                process,
                args=[file],
                kwds=kwargs
            ) for file in files]
            for file, (sample, trees[sample],
                       taxids[sample], abundances[sample],
                       accs[sample], scores[sample]
                       ) in zip(files, [r.get() for r in async_results]):
                samples.append(sample)
    else:  # sequential processing of each sample
        for file in files:
            (sample, trees[sample],
             taxids[sample], abundances[sample],
             accs[sample], scores[sample]) = process(file, **kwargs)
            samples.append(sample)
    raw_samples.extend(samples)  # Store raw sample names
    #
    # Cross analysis of samples in parallel by taxlevel
    #
    # Avoid if just a single report file of explicitly stated by flag
    if len(files) > 1 and not avoidcross:
        print('\033[90mPlease, wait. ' +
              'Performing cross analysis in parallel...\033[0m\n')
        kwargs.update({'trees': trees, 'taxids': taxids,
                       'abundances': abundances, 'files': files,
                       'control': control})
        if platform.system() and not sequential:  # Only for known platforms
            mpctx = mp.get_context('spawn')  # Important for OSX&Win
            with mpctx.Pool(processes=min(os.cpu_count(), len(
                    Rank.selected_ranks))) as pool:
                async_results = [pool.apply_async(
                    process_rank,
                    args=[level],
                    kwds=kwargs
                ) for level in Rank.selected_ranks]
                for level, (filenames, abunds, accumulators, score) in zip(
                        Rank.selected_ranks,
                        [r.get() for r in async_results]):
                    samples.extend(filenames)
                    abundances.update(abunds)
                    accs.update(accumulators)
                    scores.update(score)
        else:  # sequential processing of each selected rank
            for level in Rank.selected_ranks:
                (filenames, abunds,
                 accumulators, score) = process_rank(level, **kwargs)
                samples.extend(filenames)
                abundances.update(abunds)
                accs.update(accumulators)
                scores.update(score)
    #
    # Generate Krona plot with all the results via Krona 2.0 XML spec
    #
    print('\033[90mBuilding the taxonomy multiple tree...\033[0m', end='')
    sys.stdout.flush()
    krona: KronaTree = KronaTree(samples,
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
    polytree: MultiTree = MultiTree(samples=samples)
    polytree.grow(taxonomy=ncbi,
                  abundances=abundances,
                  accs=accs,
                  scores=scores)
    print('\033[92m OK! \033[0m')
    print('\033[90mGenerating final Krona plot...\033[0m', end='')
    sys.stdout.flush()
    polytree.toxml(taxonomy=ncbi, krona=krona)
    krona.tohtml(htmlfile, pretty=False)
    print('\033[92m OK! \033[0m')
    if _use_pandas:
        print(f'\033[90mGenerating Excel {str(excel).lower()}'
              f' summary file...\033[0m', end='')
        sys.stdout.flush()
        xlsxwriter = pd.ExcelWriter(htmlfile.split('.html')[0] + '.xlsx')
        if excel is excel.FULL:
            list_rows: List = []
            polytree.to_items(taxonomy=ncbi, items=list_rows)
            # Generate the pandas DataFrame from items and export to Excel
            iterable_1 = [samples, [COUNT, UNASSIGNED, SCORE]]
            cols1 = pd.MultiIndex.from_product(iterable_1,
                                               names=['Samples', 'Stats'])
            iterable_2 = [['Details'], ['Rank', 'Name']]
            cols2 = pd.MultiIndex.from_product(iterable_2)
            columns = cols1.append(cols2)
            df: pd.DataFrame = pd.DataFrame.from_items(list_rows,
                                                       orient='index',
                                                       columns=columns)
            df.index.names = ['TaxId']
            df.to_excel(xlsxwriter, sheet_name=str(excel))
        elif excel is excel.CMPLXCRUNCHER:
            target_ranks: List = [Rank.NO_RANK]
            if control:
                target_ranks = [Rank.SPECIES, Rank.GENUS,  # Ranks of interest
                                Rank.FAMILY, Rank.ORDER]   # for cmplxcruncher
            for rank in target_ranks:  # Once for no rank dependency (NO_RANK)
                indexes: List[int]
                sheet_name: str
                columns: List[str]
                if control:
                    indexes = [i for i in range(len(raw_samples), len(samples))
                               if (samples[i].startswith(STR_CONTROL)
                                   and rank.name.lower() in samples[i])]
                    sheet_name = f'{STR_CONTROL}_{rank.name.lower()}'
                    columns = [samples[i].split('_')[2] for i in indexes]
                else:  # No rank dependency
                    indexes = list(range(len(raw_samples)))
                    sheet_name = f'raw_samples_{rank.name.lower()}'
                    columns = [samples[i].split('_')[0] for i in indexes]
                list_rows: List = []
                polytree.to_items(taxonomy=ncbi, items=list_rows,
                                  sample_indexes=indexes)
                df: pd.DataFrame = pd.DataFrame.from_items(list_rows,
                                                           orient='index',
                                                           columns=columns)
                df.index.names = ['TaxId']
                df.to_excel(xlsxwriter, sheet_name=sheet_name)
        else:
            raise Exception(
                f'\n\033[91mERROR!\033[0m Unknown Excel option "{excel}"')
        xlsxwriter.save()
        print('\033[92m OK! \033[0m')


if __name__ == '__main__':
    main()
