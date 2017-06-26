"""
Comparative metagenomics library.

A front-end interface is provided by the companion recentrifuge module,
which will import this library.

Modules included:
    config: provides constants and other package-wide stuff
    core:

"""

import sys

__all__ = ['config', 'core', 'centrifuge', 'krona']
__author__ = 'Jose Manuel Marti'

# python
MAJOR, MINOR, *_ = sys.version_info
PYTHON_REL = (MAJOR == 3 and MINOR >= 6)
if not PYTHON_REL:
    raise ImportError('recentrifuge requires Python 3.6 or later')
