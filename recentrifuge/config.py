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
Scores = NewType('Scores', Dict[TaxId, Score])
# pylint: enable=invalid-name

# Predefined internal constants
PATH: Filename = Filename('.')
TAXDUMP_PATH: Filename = Filename('./taxdump')
NODES_FILE: Filename = Filename('nodes.dmp')
NAMES_FILE: Filename = Filename('names.dmp')
PLASMID_FILE: Filename = Filename('plasmid.names.txt')
ZIPFILE: Filename = Filename('taxdmp.zip')
JSLIB: Filename = Filename('krona.js')
HTML_SUFFIX: Filename = Filename('.rcf.html')
STR_CONTROL: str = 'CTRL'
STR_EXCLUSIVE: str = 'EXCLUSIVE'
STR_SHARED: str = 'SHARED'
STR_CONTROL_SHARED: str = 'CONTROL_SHARED'
STR_SUMMARY: str = 'SUMMARY'
DEFMINTAXA: int = 10  # min taxa to avoid collapsing one level into the parent
UNCLASSIFIED: TaxId = TaxId('0')
ROOT: TaxId = TaxId('1')
CELLULAR_ORGANISMS: TaxId = TaxId('131567')
NO_SCORE: Score = Score(None)  # Score given to taxa with no score
EPS: float = 1e-8
ADV_CTRL_MIN_SAMPLES: int = 2  # Min num of non-ctrl samples to enable advCtrl
SEVR_CONTM_MIN_RELFREQ: float = 0.01  # Min rel frequency of severe contaminant
MILD_CONTM_MIN_RELFREQ: float = 0.001  # Min rel frequency of mild contaminant


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


class Err(Enum):
    """Enumeration with error codes"""
    NO_ERROR = 0  # No error
    VOID_CTRL = 1  # A control sample is empty
    VOID_SAMPLE = 2  # A (non-control) sample is empty


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
    elif num > 1e+3:
        value = num / 1e+3
        unit = 'knt'
    else:
        value = num
        unit = 'nt'
        return f'{value} {unit}'
    return f'{value:.2f} {unit}'
