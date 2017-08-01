#!/usr/bin/env python3
"""
Get needed taxdump files from NCBI.
"""
# pylint: disable=no-name-in-module, not-an-iterable
import argparse
from ftplib import FTP
import os
import sys
from zipfile import ZipFile

from recentrifuge.config import NODES_FILE, NAMES_FILE, ZIPFILE, TAXDUMP_PATH

__version__ = '0.0.1'
__author__ = 'Jose Manuel Marti'
__date__ = 'Ago 2017'


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
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default=TAXDUMP_PATH,
        help=('path for the nodes information files (nodes.dmp and names.dmp' +
              ' from NCBI')
    )

    # Parse arguments
    args = parser.parse_args()

    # Program header
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Load NCBI nodes, names and build children
    print(f'\033[90mDownloading {ZIPFILE} from NCBI...\033[0m', end='')
    sys.stdout.flush()
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('/pub/taxonomy/')
    ftp.retrbinary('RETR ' + ZIPFILE, open(ZIPFILE, 'wb').write)
    ftp.quit()
    print('\033[92m OK! \033[0m\n')

    filezip = ZipFile(ZIPFILE)
    for filename in [NODES_FILE, NAMES_FILE]:
        print(f'\033[90mExtracting "{filename}"...', end='')
        try:
            path = filezip.extract(filename, path=args.nodespath)
        except KeyError:
            print('\n\033[91mERROR!\033[0m')
        else:
            print('\033[92m OK! \033[0m')

if __name__ == '__main__':
    main()
