"""
Comparative metagenomics library.

A front-end interface is provided by the companion recentrifuge module,
which will import this library.

Modules included:
    config: provides constants and other package-wide stuff
    core:

"""

import sys
from Bio import SeqIO
from . import LmatIO  # LMAT support
from . import CentrifugeIO  # Centrifuge support

__all__ = ['config', 'core', 'centrifuge', 'krona']
__author__ = 'Jose Manuel Marti'

# python
MAJOR, MINOR, *_ = sys.version_info
PYTHON_REL = (MAJOR == 3 and MINOR >= 6)
if not PYTHON_REL:
    raise ImportError('recentrifuge requires Python 3.6 or later')

# Addition of new custom formats to SeqIO
SeqIO._FormatToIterator["lmat"] = LmatIO.LmatOutIterator
SeqIO._FormatToIterator["centrifuge"] = CentrifugeIO.cfg_out_iterator
