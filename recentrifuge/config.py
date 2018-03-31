"""
This module provides constants and other package-wide stuff.

"""
from enum import Enum
from statistics import mean
from typing import Dict, List, Counter, NewType, Union, NamedTuple

from recentrifuge.shared_counter import SharedCounter

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
UnionCounter = Union[Counter[TaxId], SharedCounter]
UnionScores = Union[Scores, SharedCounter]
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
EPS: float = 1e-14
ROBUST_MIN_SAMPLES: int = 1  # Min num of non-ctrl samples to enable advCtrl
ROBUST_XOVER_OUTLIER = 5  # Cutoff for crossover outlier test, typ in [2, 3]
ROBUST_XOVER_ORD_MAG = 3  # Relfreq order of magnitude dif in crossover test
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


class NT(int):
    """Class representing any chain of nucleotides"""

    def __str__(self) -> str:
        """Format nucleotides number with SI prefixes and units"""
        value: float
        unit: str
        if self > 1e+12:
            value = self / 1e+12
            unit = 'Tnt'
        elif self > 1e+9:
            value = self / 1e+9
            unit = 'Gnt'
        elif self > 1e+6:
            value = self / 1e+6
            unit = 'Mnt'
        elif self > 1e+3:
            value = self / 1e+3
            unit = 'knt'
        else:
            value = self
            unit = 'nt'
            return f'{value:d} {unit}'
        return f'{value:.2f} {unit}'


# pylint: disable=too-few-public-methods
class SeqsStats(NamedTuple):
    """Sequence statistics"""
    read: int = 0
    unclas: int = 0
    clas: int = 0
    filt: int = 0


class ScoreStats(NamedTuple):
    """Score statistics"""
    maxi: Score = NO_SCORE
    mean: Score = NO_SCORE
    mini: Score = NO_SCORE


class LengthStats(NamedTuple):
    """Leg statistics"""
    maxi: NT = NT(0)
    mean: NT = NT(0)
    mini: NT = NT(0)
# pylint: enable=too-few-public-methods


class SampleStats(object):
    """Sample statistics"""

    def __init__(self, is_ctrl: bool = False,
                 minscore: Score = None, nt_read: int = 0,
                 seq_read: int = 0, seq_filt: int = 0,
                 seq_clas: int = None, seq_unclas: int = None,
                 scores: Dict[TaxId, List[Score]] = None,
                 lens: Dict[TaxId, List[int]] = None) -> None:
        """Initialize some data and setup data structures"""
        self.is_ctrl: bool = is_ctrl
        self.minscore: Score = minscore
        self.nt_read: NT = NT(nt_read)
        self.seq: SeqsStats
        if seq_clas:
            self.seq = SeqsStats(
                read=seq_read, unclas=seq_read - seq_clas,
                clas=seq_clas, filt=seq_filt)
        else:
            self.seq = SeqsStats(
                read=seq_read, unclas=seq_unclas,
                clas=seq_read - seq_unclas, filt=seq_filt)
        if scores:
            self.sco: ScoreStats = ScoreStats(
                mini=Score(min([min(s) for s in scores.values()])),
                mean=Score(mean([mean(s) for s in scores.values()])),
                maxi=Score(max([max(s) for s in scores.values()])))
        else:
            self.sco = ScoreStats()
        if lens:
            self.len: LengthStats = LengthStats(
                mini=NT(min([min(l) for l in lens.values()])),
                mean=NT(mean([mean(l) for l in lens.values()])),
                maxi=NT(max([max(l) for l in lens.values()])))
        else:
            self.len = LengthStats()
        self.num_taxa: int = len(scores)

    def to_dict(self) -> Dict[str, Union[int, Score]]:
        """
        Create a dict with the data of the object (used to feed a DataFrame)
        """
        return {'Seqs. read': self.seq.read, 'Seqs. unclass.': self.seq.unclas,
                'Seqs. class.': self.seq.clas, 'Seqs. filtered': self.seq.filt,
                'Score min': self.sco.mini, 'Score mean': self.sco.mean,
                'Score max': self.sco.maxi, 'Length min': self.len.mini,
                'Length mean': self.len.mean, 'Length max': self.len.maxi,
                'Total nt read': self.nt_read, 'Taxa assigned': self.num_taxa,
                'Score limit': self.minscore}

    def to_krona(self) -> Dict[str, str]:
        """
        Create a dict with the data of the object (used to feed a Krona plot)
        """
        return {'isctr': str(self.is_ctrl),
                'sread': str(self.seq.read),
                'sclas': str(self.seq.clas),
                'sfilt': str(self.seq.filt),
                'scmin': str(self.sco.mini),
                'scavg': str(self.sco.mean),
                'scmax': str(self.sco.maxi),
                'lnmin': str(self.len.mini),
                'lnavg': str(self.len.mean),
                'lnmax': str(self.len.maxi),
                'taxas': str(self.num_taxa),
                'sclim': str(self.minscore),
                'totnt': str(self.nt_read)}

    def get_unclas_ratio(self) -> float:
        """Get ratio of unclassified sequences"""
        return self.seq.unclas / self.seq.read

    def get_reject_ratio(self) -> float:
        """Get ratio of rejected sequences by filtering"""
        return 1 - self.seq.filt / self.seq.clas

def ansi(num: int):
    """Return function that escapes text with ANSI color n."""
    return lambda txt: f'\033[{num}m{txt}\033[0m'


# pylint: disable=invalid-name
gray, red, green, yellow, blue, magenta, cyan, white = map(ansi, range(90, 98))
# pylint: enable=invalid-name
