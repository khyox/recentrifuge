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
# pylint: enable=invalid-name

# Predefined internal constants
PATH = '.'
NODESFILE = 'nodes.dmp'
NAMESFILE = 'names.dmp'
LINEAGEXT = '.lineage'
HTML_SUFFIX = '.rcf.html'
DEFMINTAXA = 10
ROOT = TaxId('1')
