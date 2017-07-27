#!/usr/bin/env python3
"""
Extract reads following Centrifuge/Kraken output.
"""
# pylint: disable=no-name-in-module, not-an-iterable
import argparse
import os
import sys
from collections import Counter
from typing import List, Set

from Bio import SeqIO, SeqRecord

from recentrifuge.config import Filename, TaxId
from recentrifuge.config import NODESFILE, NAMESFILE
from recentrifuge.core import Taxonomy, TaxLevels, TaxTree, Rank, Ranks

__version__ = '0.0.1'
__author__ = 'Jose Manuel Marti'
__date__ = 'Jul 2017'


def main():
    """Main entry point to script."""
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Extract reads following Centrifuge/Kraken output',
        epilog=f'%(prog)s  - {__author__} - {__date__}'
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s release {__version__} ({__date__})'
    )
    parser.add_argument(
        '-f', '--file',
        action='store',
        metavar='FILE',
        help='Centrifuge output file.'
    )
    parser.add_argument(
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default='./',
        help=('path for the nodes information files (nodes.dmp and names.dmp' +
              ' from NCBI')
    )
    parser.add_argument(
        '-i', '--include',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to include a taxon and all underneath ' +
              '(multiple -i is available to include several taxid). ' +
              'By default all the taxa is considered for inclusion.')
    )
    parser.add_argument(
        '-x', '--exclude',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to exclude a taxon and all underneath ' +
              '(multiple -x is available to exclude several taxid)')
    )
    parser.add_argument(
        '-1', '--mate1',
        action='store',
        metavar='FILE',
        default=None,
        help='Paired-ends FASTQ file for mate 1s '
             '(filename usually includes _1)'
    )
    parser.add_argument(
        '-2', '--mate2',
        action='store',
        metavar='FILE',
        default=None,
        help='Paired-ends FASTQ file for mate 2s '
             '(filename usually includes _2)'
    )

    # Parse arguments
    args = parser.parse_args()
    output_file = args.file
    nodesfile: Filename = Filename(os.path.join(args.nodespath, NODESFILE))
    namesfile: Filename = Filename(os.path.join(args.nodespath, NAMESFILE))
    excluding: Set[TaxId] = set(args.exclude)
    including: Set[TaxId] = set(args.include)
    fastq_1: Filename = args.mate1
    fastq_2: Filename = args.mate2

    # Program header and chdir
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Load NCBI nodes, names and build children
    ncbi: Taxonomy = Taxonomy(nodesfile, namesfile,
                              False, excluding, including)

    # Build taxonomy tree
    print('\033[90mBuilding taxonomy tree...\033[0m', end='')
    sys.stdout.flush()
    tree = TaxTree()
    tree.grow(taxonomy=ncbi)
    print('\033[92m OK! \033[0m')

    # Get the taxa
    print('\033[90mFiltering taxa...\033[0m', end='')
    sys.stdout.flush()
    ranks: Ranks = Ranks({})
    tree.get_taxa(ranks=ranks,
                  include=including,
                  exclude=excluding)
    print('\033[92m OK! \033[0m')
    taxids: Set[TaxId] = set(ranks)
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    num_taxlevels = Counter({rank: len(taxlevels[rank]) for rank in taxlevels})
    num_taxlevels = +num_taxlevels

    # Statistics about including taxa
    print(f'  {len(ranks)}\033[90m taxid selected in \033[0m', end='')
    print(f'{len(num_taxlevels)}\033[90m different taxonomical levels:\033[0m')
    for rank in num_taxlevels:
        print(f'  Number of different {rank}: {num_taxlevels[rank]}')

    # Get the records
    records: List[SeqRecord] = []
    num_seqs: int = 0
    print(f'\033[90mLoading output file {output_file}...\033[0m', end='')
    sys.stdout.flush()
    try:
        with open(output_file, 'rU') as file:
            file.readline()  # discard header
            for num_seqs, record in enumerate(SeqIO.parse(file, 'centrifuge')):
                tid: TaxId = record.annotations['taxID']
                if tid not in taxids:
                    continue
                # score: Score = record.annotations['score']
                records.append(record)
    except FileNotFoundError:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        output_file + '"')
    print('\033[92m OK! \033[0m')

    # Basic records statistics
    print(f'  \033[90mMatching reads: \033[0m{len(records):_d} \033[90m\t'
          f'(\033[0m{len(records)/num_seqs:.4%}\033[90m of sample)')

    # FASTQ sequence dealing
    print(f'\033[90mLoading FASTQ files {fastq_1} and {fastq_2}...\n'
          f'Mseqs: \033[0m', end='')
    sys.stdout.flush()
    records_ids: List[SeqRecord] = [record.id for record in records]
    seqs1: List[SeqRecord] = []
    seqs2: List[SeqRecord] = []
    try:
        with open(fastq_1, 'rU') as file1, open(fastq_2, 'rU') as file2:
            for i, (rec1, rec2) in enumerate(zip(SeqIO.parse(file1, 'fastq'),
                                                 SeqIO.parse(file2, 'fastq'))):
                if not i % 1000000:
                    print(f'{i//1000000:_d}', end='')
                    sys.stdout.flush()
                elif not i % 100000:
                    print('.', end='')
                    sys.stdout.flush()
                try:
                    records_ids.remove(rec1.id)
                except ValueError:
                    pass
                else:
                    seqs1.append(rec1)
                    seqs2.append(rec2)
    except FileNotFoundError:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        output_file + '"')
    print('\033[92m OK! \033[0m')
    filename1 = f'{fastq_1}.{".".join(taxids)}.fastq'
    filename2 = f'{fastq_2}.{".".join(taxids)}.fastq'
    SeqIO.write(seqs1, filename1, 'fastq')
    print(f'\033[90mWrote \033[0m{len(seqs1)}\033[90m reads in {filename1}')
    SeqIO.write(seqs2, filename2, 'fastq')
    print(f'\033[90mWrote \033[0m{len(seqs2)}\033[90m reads in {filename2}')


if __name__ == '__main__':
    main()
