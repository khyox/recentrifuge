"""
This module provides constants and other package-wide stuff.

"""
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
NODESFILE = 'nodes.dmp'
NAMESFILE = 'names.dmp'
HTML_SUFFIX = '.rcf.html'
DEFMINTAXA = 10  # minimum taxa to avoid collapsing one level to the parent one
UNCLASSIFIED = TaxId('0')
ROOT = TaxId('1')
CELLULAR_ORGANISMS = TaxId('131567')
NO_SCORE = Score(0)  # score given to taxa with no score
