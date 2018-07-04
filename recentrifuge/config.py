"""
This module provides constants and other package-wide stuff.

"""
from enum import Enum
from typing import Dict, Counter, NewType, Union

from recentrifuge.shared_counter import SharedCounter

# Type annotations
# pylint: disable=invalid-name
Filename = NewType('Filename', str)
Sample = NewType('Sample', str)
Id = str
Children = NewType('Children', Dict[Id, Dict[Id, int]])
Names = NewType('Names', Dict[Id, str])
Parents = NewType('Parents', Dict[Id, Id])
Score = NewType('Score', float)
Scores = NewType('Scores', Dict[Id, Score])
UnionCounter = Union[Counter[Id], SharedCounter]
UnionScores = Union[Scores, SharedCounter]
# pylint: enable=invalid-name

# Predefined internal constants
PATH: Filename = Filename('.')
GODUMP_PATH: Filename = Filename('./godump')
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
DEFMINGENE: int = 2  # min genes to avoid collapsing GO genes hierarchy levels
ROOT: Id = Id('1')
CELLULAR_ORGANISMS: Id = Id('131567')
EPS: float = 1e-14
ROBUST_MIN_SAMPLES: int = 1  # Min num of non-ctrl samples to enable advCtrl
ROBUST_XOVER_OUTLIER = 5  # Cutoff for crossover outlier test, typ in [2, 3]
ROBUST_XOVER_ORD_MAG = 3  # Relfreq order of magnitude dif in crossover test
SEVR_CONTM_MIN_RELFREQ: float = 0.01  # Min rel frequency of severe contaminant
MILD_CONTM_MIN_RELFREQ: float = 0.001  # Min rel frequency of mild contaminant
GO_ROOT: Id = Id('GO:ROOT')
NO_SCORE = Score(float("-inf"))


class Classifier(Enum):
    """Enumeration with supported taxonomic classifiers."""
    CENTRIFUGE = 0
    LMAT = 1
    CLARK = 2
    KRAKEN = 3

    def __str__(self):
        return f'{str(self.name)}'


class Chart(Enum):
    """Enumeration with Krona chart options."""
    TAXOMIC = 0  # Recentrifuge plot
    GENOMIC = 1  # Regentrifuge plot

    def __str__(self):
        return f'{str(self.name)}'


class Scoring(Enum):
    """Enumeration with scoring options."""
    SHEL = 0  # Single Hit Equivalent Length (default)
    LENGTH = 1  # The length of a read or the combined length of mate pairs
    LOGLENGTH = 2  # Log10 of the length
    NORMA = 3  # It is the normalized score SHEL / LENGTH
    LMAT = 4  # Default for LMAT, not available for others
    CLARK = 5  # Based on CLARK's confidence score, not available for others

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


class Unscore(Enum):
    """Enumeration with special score cases (provisional)."""
    NO_SCORE = 0  # Given to taxa with no score

    def __str__(self):
        return f'{str(self.name)}'

# NO_SCORE = Unscore.NO_SCORE
