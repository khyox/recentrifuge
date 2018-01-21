"""
Comparative metagenomics library.

A front-end interface is provided by the companion recentrifuge module,
which will import this library.

"""

import sys
from Bio import SeqIO
from . import lmat_io  # LMAT support
from . import centrifuge_io  # Centrifuge support

__all__ = ['config', 'shared_counter', 'taxonomy', 'trees', 'rank',
           'core', 'centrifuge', 'lmat', 'krona']
__author__ = 'Jose Manuel Marti'

# python
MAJOR, MINOR, *_ = sys.version_info
PYTHON_REL = (MAJOR == 3 and MINOR >= 6)
if not PYTHON_REL:
    raise ImportError('Recentrifuge requires Python 3.6 or later')

# Addition of new custom formats to SeqIO
# pylint: disable=protected-access
SeqIO._FormatToIterator["lmat"] = lmat_io.lmat_out_iterator
SeqIO._FormatToIterator["centrifuge"] = centrifuge_io.cfg_out_iterator
# pylint: enable=protected-access
