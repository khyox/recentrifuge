"""
This module provides constants and other package-wide stuff.

"""
from enum import Enum
from typing import Dict, Counter, NewType, Union

from recentrifuge.shared_counter import SharedCounter

# Licence notice
LICENSE = ('''

    Copyright (C) 2017, 2018, Jose Manuel Martí Martínez
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
''')

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
XLSX_SUFFIX: Filename = Filename('.rcf.xlsx')
MOCK_XLSX: Filename = Filename('mock.xlsx')
MOCK_SHEET: Filename = Filename('mock')
KEY_SHEET: Filename = Filename('key')
REXTRACT_TEST_SAMPLE: Filename = Filename('smpl3.out')
REXTRACT_TEST_FASTQ: Filename = Filename('smpl3.fastq')
REXTRACT_TEST_TAXID: Id = Id('2208')
TEST_INPUT_DIR: Filename = Filename('test/')
TEST_OUTPUT_DIR: Filename = Filename('rcf_test/')
STATS_SHEET_NAME: str = '_sample_stats'
STR_CONTROL: str = 'CTRL'
STR_EXCLUSIVE: str = 'EXCLUSIVE'
STR_SHARED: str = 'SHARED'
STR_CONTROL_SHARED: str = 'CONTROL_SHARED'
STR_SUMMARY: str = 'SUMMARY'
DEFMINGENE: int = 2  # min genes to avoid collapsing GO genes hierarchy levels
ROOT: Id = Id('1')
CELLULAR_ORGANISMS: Id = Id('131567')
GO_ROOT: Id = Id('GO:ROOT')
NO_SCORE = Score(float("-inf"))


class Classifier(Enum):
    """Enumeration with supported taxonomic classifiers."""
    CENTRIFUGE = 0
    LMAT = 1
    CLARK = 2
    KRAKEN = 3
    GENERIC = 4

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
    NORMA = 3  # This is the normalized score SHEL / LENGTH
    LMAT = 4  # Default for LMAT, not available for others
    CLARK_C = 5  # Based on CLARK's confidence score, not available for others
    CLARK_G = 6  # Based on CLARK's gamma score, not available for others
    KRAKEN = 7  # % of k-mers that mapped to the taxid labeled in kraken
    GENERIC = 8  # Score from a generic classifier

    def __str__(self):
        return f'{self.name}'


class Extra(Enum):
    """Enumeration with options for the extra output."""
    FULL = 0  # Provide excel with detailed results including score (default)
    CMPLXCRUNCHER = 1  # Results in excel with cmplxcruncher format
    CSV = 2  # Results as a collection of CSV files
    TSV = 3  # Results as a collection of TSV files

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
