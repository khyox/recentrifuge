"""
Rank class for representing taxonomic levels or ranks.

"""

from enum import Enum
from typing import List, Iterator, NewType, Dict, Set

from recentrifuge.config import Id

# Type annotations
# pylint: disable=invalid-name
# Ranks and Levels are devised to be one the inverse of the other
Ranks = NewType('Ranks', Dict[Id, 'Rank'])  # Rank of each Id
TaxLevels = NewType('TaxLevels', Dict['Rank', Set[Id]])  # TaxIds for rank
# pylint: enable=invalid-name


class classproperty(object):  # pylint: disable=invalid-name
    """Decorator to emulate a class property."""

    # pylint: disable=too-few-public-methods
    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)


class UnsupportedTaxLevelError(Exception):
    """Raised if a unsupported tax level is found."""


class Rank(Enum):
    """Enumeration with ranks (taxonomical levels).

    The members are initialized using a customized __new__ method so
    that there are 3 different mechanisms to assign value to members:
      1) Empty will imply decrease autonumbering starting at 99.
      2) Direct assignation of an integer value is allowed.
      3) A string with the name of a previously listed member will
        assign that same value to the current member.
    With this, different comparisons between ranks will be available.

    """
    # pylint: disable=invalid-name
    ROOT = ()
    DOMAIN = ()
    D = 'DOMAIN'
    SUPERKINGDOM = ()
    KINGDOM = ()
    K = 'KINGDOM'
    SUBKINGDOM = ()
    SUPERPHYLUM = ()
    PHYLUM = ()
    P = 'PHYLUM'
    SUBPHYLUM = ()
    SUPERCLASS = ()
    CLASS = ()
    C = 'CLASS'
    SUBCLASS = ()
    INFRACLASS = ()
    COHORT = ()
    SUBCOHORT = ()
    SUPERORDER = ()
    ORDER = ()
    O = 'ORDER'
    SUBORDER = ()
    INFRAORDER = ()
    PARVORDER = ()
    SUPERFAMILY = ()
    FAMILY = ()
    F = 'FAMILY'
    SUBFAMILY = ()
    TRIBE = ()
    SUBTRIBE = ()
    GENUS = ()
    G = 'GENUS'
    SUBGENUS = ()
    SECTION = ()
    SUBSECTION = ()
    SERIES = ()
    SPECIES_GROUP = ()
    SPECIES_SUBGROUP = ()
    SPECIES = ()
    S = 'SPECIES'
    SUBSPECIES = ()
    VARIETAS = ()
    FORMA = ()
    GO0 = ()
    GO1 = ()
    GO2 = ()
    GO3 = ()
    GO4 = ()
    GO5 = ()
    GO6 = ()
    GO7 = ()
    GO8 = ()
    GO9 = ()
    UNCLASSIFIED = 0
    U = 'UNCLASSIFIED'
    NO_RANK = -1
    # pylint: enable=invalid-name

    @classproperty
    def selected_ranks(cls):  # pylint: disable=no-self-argument
        """Ranks selected for deep analysis and comparisons"""
        _selected_taxlevels: List['Rank'] = [cls.S, cls.G, cls.F, cls.O,
                                             cls.C, cls.P, cls.K, cls.D]
#        _selected_taxlevels: List['Rank'] = [cls.S, cls.G, cls.F, cls.O]
        return _selected_taxlevels

    @classproperty
    def genomic_ranks(cls):  # pylint: disable=no-self-argument
        """GO ranks selected for deep analysis and comparisons"""
        _selected_golevels: List['Rank'] = [cls.GO9, cls.GO8, cls.GO7, cls.GO6,
                                            cls.GO5, cls.GO4, cls.GO3, cls.GO2]
        return _selected_golevels

    @classmethod
    def centrifuge(cls, tax_level: str) -> 'Rank':
        """Transforms Centrifuge codes for taxonomical levels"""
        taxonomic_level: Rank
        if tax_level == '-':
            taxonomic_level = cls.NO_RANK
        else:
            try:
                taxonomic_level = cls[tax_level]
            except KeyError:
                raise UnsupportedTaxLevelError(
                    f'Unknown tax level {tax_level}')
        return taxonomic_level

    @classmethod
    def ranks_to_taxlevels(cls, ranks: Ranks) -> TaxLevels:
        """Generate TaxLevels (taxids of ranks) from Ranks (rank of taxids)."""
        return TaxLevels({rank: {taxid for taxid in ranks if
                                 ranks[taxid] is rank} for rank in Rank})

    @property
    def ranks_from_specific(self) -> Iterator['Rank']:
        """Generator returning selected taxlevels from specific to general."""
        for rank in list(Rank.selected_ranks):
            yield rank
            if rank is self:
                break

    @property
    def ranks_from_general(self) -> Iterator['Rank']:
        """Generator returning selected taxlevels from general to specific."""
        for rank in reversed(list(Rank.selected_ranks)):
            yield rank
            if rank is self:
                break

    def __new__(cls, init=None):
        _value: int
        if isinstance(init, int):
            _value = init
        elif isinstance(init, str):
            _value = cls[init].value
        else:
            _value = 99 - len(cls.__members__)  # Give decreasing values < 100
        obj = object.__new__(cls)
        obj._value_ = _value  # pylint: disable=protected-access
        return obj

    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)

    def __str__(self):
        return f'{self.name}'

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            if self.value > 0 and other.value > 0:
                return self.value < other.value
            return False
        return NotImplemented

    def __le__(self, other):
        if self.__class__ is other.__class__:
            if self.value > 0 and other.value > 0:
                return self.value <= other.value
            return False
        return NotImplemented

    def __gt__(self, other):
        if self.__class__ is other.__class__:
            if self.value > 0 and other.value > 0:
                return self.value > other.value
            return False
        return NotImplemented

    def __ge__(self, other):
        if self.__class__ is other.__class__:
            if self.value > 0 and other.value > 0:
                return self.value >= other.value
            return False
        return NotImplemented
