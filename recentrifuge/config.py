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
