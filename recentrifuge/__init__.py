#
#     Copyright (C) 2017, 2018, Jose Manuel Martí Martínez
#
#     With the single exception of the file krona.js, which has embedded
#     its own free software licence, you can redistribute Recentrifuge
#     and/or modify it under the terms of the GNU Affero General Public
#     License as published by the Free Software Foundation, either
#     version 3 of the License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Robust comparative analysis and contamination removal for metagenomics library

A front-end interface is provided by the companion Recentrifuge module,
which will import this library.

"""

__all__ = ['config', 'shared_counter', 'taxonomy', 'trees', 'rank', 'core',
           'centrifuge', 'lmat', 'clark', 'krona', 'kraken', 'generic',
           'ontology', 'params', 'stats', 'taxclass', 'mock',
           '__author__', '__date__', '__version__']
__author__ = 'Jose Manuel Martí'
__copyright__ = 'Copyright (C) 2017–2019 Jose Manuel Martí Martínez'
__license__ = 'GNU Affero General Public License Version 3'
__email__ = 'jse.mnl **AT** gmail.com'
__maintainer__ = 'Jose Manuel Martí'
__status__ = 'Beta'
__date__ = 'May 2019'
__version__ = '0.28.9'

import sys
from Bio import SeqIO
from . import lmat_io  # LMAT support
from . import centrifuge_io  # Centrifuge support
from . import fastq_io  # Quick FASTQ support

# python
MAJOR, MINOR, *_ = sys.version_info
PYTHON_REL = (MAJOR == 3 and MINOR >= 6)
if not PYTHON_REL:
    raise ImportError('Recentrifuge requires Python 3.6 or later')

# Addition of new custom formats to SeqIO
# pylint: disable=protected-access
SeqIO._FormatToIterator["lmat"] = lmat_io.lmat_out_iterator
SeqIO._FormatToIterator["centrifuge"] = centrifuge_io.cfg_out_iterator
SeqIO._FormatToIterator["quickfastq"] = fastq_io.quick_fastq_iterator
SeqIO._FormatToWriter["lmat"] = lmat_io.LmatOutWriter
SeqIO._FormatToWriter["quickfastq"] = fastq_io.QuickFastqWriter
# pylint: enable=protected-access
