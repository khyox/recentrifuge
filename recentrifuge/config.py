"""
This module provides constants and other package-wide stuff.

"""
from enum import Enum
from typing import Dict, NewType

# Type annotations
# pylint: disable=invalid-name
Filename = NewType('Filename', str)
Sample = NewType('Sample', str)
TaxId = NewType('TaxId', str)
Children = NewType('Children', Dict[TaxId, Dict[TaxId, int]])
Names = NewType('Names', Dict[TaxId, str])
Parents = NewType('Parents', Dict[TaxId, TaxId])
Score = NewType('Score', float)
# pylint: enable=invalid-name

# Predefined internal constants
PATH = '.'
TAXDUMP_PATH = './taxdump'
NODES_FILE = 'nodes.dmp'
NAMES_FILE = 'names.dmp'
PLASMID_FILE = 'plasmid.names.txt'
ZIPFILE = 'taxdmp.zip'
JSLIB = 'krona.js'
HTML_SUFFIX = '.rcf.html'
STR_CONTROL = 'CTRL'
STR_EXCLUSIVE = 'EXCLUSIVE'
STR_SHARED = 'SHARED'
STR_SHARED_CONTROL = 'SHARED_CONTROL'
DEFMINTAXA = 10  # minimum taxa to avoid collapsing one level to the parent one
UNCLASSIFIED = TaxId('0')
ROOT = TaxId('1')
CELLULAR_ORGANISMS = TaxId('131567')
NO_SCORE = Score(None)  # Score given to taxa with no score
EPS = 1e-8

class Scoring(Enum):
    """Enumeration with scoring options."""
    SHEL = 0  # Single Hit Equivalent Length (default)
    LENGTH = 1  # The length of a read or the combined length of mate pairs
    LOGLENGTH = 2  # Log10 of the length
    NORMA = 3  # It is the normalized score SHEL / LENGTH
    LMAT = 4  # Default for LMAT, not available for Centrifuge

    def __str__(self):
        return f'{self.name}'


class Excel(Enum):
    """Enumeration with excel output options."""
    FULL = 0  # Provide detailed results including score (default)
    CMPLXCRUNCHER = 1  # Results in cmplxcruncher format

    def __str__(self):
        return f'{str(self.name)}'


def ansi(num: int):
    """Return function that escapes text with ANSI color n."""
    return lambda txt: f'\033[{num}m{txt}\033[0m'


# pylint: disable=invalid-name
gray, red, green, yellow, blue, magenta, cyan, white = map(ansi, range(90, 98))
# pylint: enable=invalid-name


def nucleotides(num: int) -> str:
    """Format nucleotides number with SI prefixes and units"""
    value: float
    unit: str
    if num > 1e+12:
        value = num / 1e+12
        unit = 'Tnt'
    elif num > 1e+9:
        value = num / 1e+9
        unit = 'Gnt'
    elif num > 1e+6:
        value = num / 1e+6
        unit = 'Mnt'
    elif num > 1e+6:
        value = num / 1e+3
        unit = 'knt'
    else:
        value = num
        unit = 'nt'
    return f'{value:.2f} {unit}'
