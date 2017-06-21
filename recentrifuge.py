#!/usr/bin/env python3
"""
Post-process Centrifuge/Kraken report files for Krona.
"""
# pylint: disable=not-an-iterable, bad-reversed-sequence, E0203
__version__ = '0.6.3'
__author__ = 'Jose Manuel Marti'
__date__ = 'Jun 2017'

import argparse
import collections as col
import copy
import csv
import io
import multiprocessing as mp
import os
import platform
import subprocess
import sys
from enum import Enum
from typing import Iterable, Dict, Counter, Tuple, List, Set
from typing import NewType, Iterator, Union, Optional

# Type annotations
# pylint: disable=invalid-name
Filename = NewType('Filename', str)
Sample = NewType('Sample', str)
TaxId = NewType('TaxId', str)
Children = NewType('Children', Dict[TaxId, Dict[TaxId, int]])
Names = NewType('Names', Dict[TaxId, str])
Parents = NewType('Parents', Dict[TaxId, TaxId])
# Ranks and Levels are devised to be one the inverse of the other
Ranks = NewType('Ranks', Dict[TaxId, 'Rank'])  # Rank of each TaxId
TaxLevels = NewType('TaxLevels', Dict['Rank', Set[TaxId]])  # TaxIds for rank
# pylint: enable=invalid-name

# Predefined internal constants
_PATH = '.'
_NODESFILE = 'nodes.dmp'
_NAMESFILE = 'names.dmp'
_LINEAGEXT = '.lineage'
_HTML_SUFFIX = '.commet.html'
_DEFMINTAXA = 10
_ROOT = TaxId('1')

# Define encoding dialect for TSV files expected by Krona
csv.register_dialect('krona', 'unix', delimiter='\t', quoting=csv.QUOTE_NONE)


class classproperty(object):  # pylint: disable=invalid-name
    """Decorator to emulate a class property."""

    # pylint: disable=too-few-public-methods
    def __init__(self, fget):
        self.fget = fget

    def __get__(self, owner_self, owner_cls):
        return self.fget(owner_cls)


class UnsupportedTaxLevelError(Exception):
    """Raised if a unsupported tx level is found."""


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
    SPECIES_GROUP = ()
    SPECIES_SUBGROUP = ()
    SPECIES = ()
    S = 'SPECIES'
    SUBSPECIES = ()
    VARIETAS = ()
    FORMA = ()
    UNCLASSIFIED = 0
    U = 'UNCLASSIFIED'
    NO_RANK = -1

    # pylint: enable=invalid-name

    @classproperty
    def selected_ranks(cls):  # pylint: disable=no-self-argument
        """Ranks selected for deep analysis and comparisons"""
        _selected_taxlevels: List['Rank'] = [cls.S, cls.G, cls.F, cls.O,
                                             cls.C, cls.P, cls.K, cls.D]
        return _selected_taxlevels

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
        for rank in Rank.selected_ranks:
            yield rank
            if rank is self:
                break

    @property
    def ranks_from_general(self) -> Iterator['Rank']:
        """Generator returning selected taxlevels from general to specific."""
        for rank in reversed(Rank.selected_ranks):
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
        obj._value_ = _value
        return obj

    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)

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


class Taxonomy:
    """Taxonomy related data and methods."""

    def __init__(self,
                 nodes_file: str,
                 names_file: str,
                 collapse: bool,
                 excluding: Set[TaxId],
                 including: Set[TaxId],
                 ) -> None:

        # Type data declaration and initialization
        self.parents: Parents = Parents({})
        self.ranks: Ranks = Ranks({})
        self.names: Names = Names({})
        self.children: Children = Children({})
        self.collapse: bool = collapse

        # Initialization methods
        self.read_nodes(nodes_file)
        self.read_names(names_file)
        self.build_children()

        # Show explicitly included and excluded taxa
        if including:
            print('List of taxa (and below) to be explicitly included:')
            print('\t\tTaxId\tScientific Name')
            for taxid in including:
                print(f'\t\t{taxid}\t{self.names[taxid]}')
        else:
            # To excluding to operate not on single taxa but on subtrees
            including = {_ROOT}
        self.including: Set[TaxId] = including
        if excluding:
            print('List of taxa (and below) to be excluded:')
            print('\t\tTaxId\tScientific Name')
            for taxid in excluding:
                print(f'\t\t{taxid}\t{self.names[taxid]}')
        self.excluding: Set[TaxId] = excluding

    def read_nodes(self, nodes_file: str) -> None:
        """Build dicts of parent and rank for a given taxid (key)"""
        print('\033[90mLoading NCBI nodes...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(nodes_file, 'r') as file:
                for line in file:
                    _tid, _parent, _rank, *_ = line.split('\t|\t')
                    tid = TaxId(_tid)
                    parent = TaxId(_parent)
                    self.parents[tid] = parent
                    rank: Rank
                    try:
                        rank = Rank[_rank.upper().replace(" ", "_")]
                    except KeyError:
                        raise UnsupportedTaxLevelError(
                            f'Unknown tax level {_rank}')
                    self.ranks[tid] = rank

        except OSError:
            print(f'\n\033[91mERROR!\033[0m Cannot read {nodes_file}')
            raise
        else:
            print('\033[92m OK! \033[0m')

    def read_names(self, names_file: str) -> None:
        """Build dict with name for a given taxid (key)."""
        print('\033[90mLoading NCBI names...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(names_file, 'r') as file:
                for line in file:
                    if 'scientific name' in line:
                        tid, scientific_name, *_ = line.split('\t|\t')
                        self.names[TaxId(tid)] = scientific_name
        except:
            raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                            names_file + '"')
        else:
            print('\033[92m OK! \033[0m')

    def build_children(self) -> None:
        """Build dict of children for a given parent taxid (key)."""
        print('\033[90mBuilding dict of parent to children taxa...\033[0m',
              end='')
        sys.stdout.flush()
        self.children: Children = Children({})
        for tid in self.parents:
            if self.parents[tid] not in self.children:
                self.children[self.parents[tid]] = {}
            self.children[self.parents[tid]][tid] = 0
        print('\033[92m OK! \033[0m')

    def get_rank(self, taxid: TaxId) -> Rank:
        """Retrieve the rank for a TaxId."""
        return self.ranks.get(taxid, Rank.UNCLASSIFIED)


class TaxTree(dict):
    """Nodes of a taxonomical tree"""

    def __init__(self, *args,
                 counts: int = 0,
                 rank: Rank = Rank.UNCLASSIFIED
                 ) -> None:
        super().__init__(args)
        self.counts: int = counts
        self.taxlevel: Rank = rank

    def __str__(self, num_min: int = 1) -> None:
        """Recursively print populated nodes of the taxonomy tree"""
        for tid in self:
            if self[tid].counts >= num_min:
                print(f'{tid}[{self[tid].counts}]', end='')
            if self[tid]:  # Has descendants
                print('', end='->(')
                self[tid].__str__(num_min=num_min)
            else:
                print('', end=',')
        print(')', end='')

    def grow(self,
             taxonomy: Taxonomy,
             abundance: Counter[TaxId],
             taxid: TaxId = _ROOT,
             _path: List[TaxId] = None,
             ) -> None:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            abundance: counter for taxids with their abundances.
            taxid: A taxid that is the root one by default
            _path: list used by the recursive algorithm to avoid loops

        Returns: None

        """
        if not _path:
            _path = []
        if taxid not in _path:  # Avoid loops for repeated taxid (like root)
            self[taxid] = TaxTree(counts=abundance.get(taxid, 0),
                                  rank=taxonomy.get_rank(taxid))
            if taxid in taxonomy.children:  # taxid has children
                for child in taxonomy.children[taxid]:
                    self[taxid].grow(taxonomy, abundance, child,
                                     _path + [taxid])

    def trace(self,
              target: TaxId,
              nodes: List[TaxId],
              ) -> bool:
        """
        Recursively get a list of nodes from self to target taxid.

        Args:
            target: TaxId of the node to trace
            nodes: Input/output list of TaxIds: at 1st entry is empty.

        Returns:
            Boolean about the success of the tracing.

        """
        target_found: bool = False
        for tid in self:
            if self[tid]:
                nodes.append(tid)
                if target in self[tid] and target != _ROOT:
                    # Avoid to append more than once the root node
                    nodes.append(target)
                    target_found = True
                    break
                if self[tid].trace(target, nodes):
                    target_found = True
                    break
                else:
                    nodes.pop()
        return target_found

    def get_lineage(self,
                    parents: Parents,
                    taxids: Iterable,
                    ) -> Tuple[str, Dict[TaxId, List[TaxId]]]:
        """
        Build dict with taxid as keys, whose values are the list of
            nodes in the tree in the path from root to such a taxid.

        Args:
            parents: dictionary of taxids parents.
            taxids: collection with the taxids to process.

        Returns:
            Log string, dict with the list of nodes for each taxid.

        """
        output = io.StringIO(newline='')
        output.write('  \033[90mGetting lineage of taxa...\033[0m')
        nodes_traced: Dict[TaxId, List[TaxId]] = {}
        for tid in taxids:
            if tid == _ROOT:
                nodes_traced[_ROOT] = [_ROOT, ]  # Root node special case
            elif tid in parents:
                nodes: List[TaxId] = []
                if self.trace(tid, nodes):  # list nodes is populated
                    nodes_traced[tid] = nodes
                else:
                    output.write('[\033[93mWARNING\033[0m: Failed tracing '
                                 f'of taxid {tid}: missing in tree]\n')
            else:
                output.write('[\033[93mWARNING\033[0m: Discarded unknown '
                             f'taxid {tid}: missing in parents]\n')
        output.write('\033[92m OK! \033[0m\n')
        return output.getvalue(), nodes_traced

    def get_taxa(self,
                 abundance: Counter[TaxId],
                 ranks: Ranks = None,
                 mindepth: int = 0,
                 maxdepth: int = 0,
                 include: Union[Tuple, Set[TaxId]] = (),
                 exclude: Union[Tuple, Set[TaxId]] = (),
                 just_level: Rank = None,
                 _in_branch: bool = False
                 ) -> None:
        """
        Recursively get the taxa between min and max depth levels.

        Args:
            abundance: Input/output Counter; at 1st entry it's empty.
            ranks: Optional I/O dict; at 1st entry should be empty.
            mindepth: 0 gets the taxa from the very beginning depth.
            maxdepth: 0 does not stop the search at any depth level.
            include: contains the root taxid of the subtrees to be
                included. If it is empty (default) all the taxa is
                included (except explicitly excluded).
            exclude: contains the root taxid of the subtrees to be
                excluded
            just_level: If set, just taxa in this taxlevel will be
                counted.
            _in_branch: is like a static variable to tell the
                recursive function that it is in a subtree that is
                a branch of a taxon in the include list.

        Returns: None

        """
        mindepth -= 1
        if maxdepth != 1:
            maxdepth -= 1
            for tid in self:
                in_branch: bool = (
                    (_in_branch or  # called with flag activated? or
                     not include or  # include by default? or
                     (tid in include))  # tid is to be included?
                    and tid not in exclude  # and not in exclude list
                )
                if (mindepth <= 0
                    and in_branch
                    and (just_level is None
                         or self[tid].taxlevel is just_level)):
                    abundance[tid] = self[tid].counts
                    if ranks is not None:
                        ranks[tid] = self[tid].taxlevel
                if self[tid]:
                    self[tid].get_taxa(abundance, ranks,
                                       mindepth, maxdepth,
                                       include, exclude,
                                       just_level, in_branch)

    def prune(self,
              mintaxa: int = 1,
              minlevel: Rank = None,
              collapse: bool = True,
              verb: bool = False,
              ) -> bool:
        """
        Recursively prune/collapse low abundant taxa of the TaxTree.

        Args:
            mintaxa: minimum taxa to avoid pruning/collapsing
                one level to the parent one.
            minlevel: if any, minimum Rank allowed in the TaxTree.
            collapse: selects if a lower level should be accumulated in
                the higher one before pruning a node (do so by default).
            verb: increase output verbosity

        Returns: True if this node is a leaf

        """
        for tid in list(self):
            if (self[tid]  # If the subtree has branches,
                and self[tid].prune(  # then prune that subtree
                    mintaxa, minlevel, collapse, verb)):
                # The pruned subtree keeps having branches (still not a leaf!)
                if verb:
                    print(f'NOT pruning branch {tid} with {self[tid].counts}')
            else:  # self[tid] is a leaf (after pruning its branches or not)
                if (self[tid].counts < mintaxa  # Not enough counts, or check
                    or (minlevel  # if minlevel is set, then check if level
                        and (self[tid].taxlevel < minlevel  # is lower or
                             or self.taxlevel <= minlevel))):  # other test
                    if collapse:
                        self.counts += self[
                            tid].counts  # Accumulate abundance in higher tax
                    if verb and self[tid].counts:
                        print(f'Pruning branch {tid} with {self[tid].counts}')
                    self.pop(tid)  # Prune leaf
                elif verb:
                    print(f'NOT pruning leaf {tid} with {self[tid].counts}')
        return bool(self)  # True if this node has branches (is not a leaf)


class SharedCounter(col.Counter):
    """Extends collection.Counter with useful ops. for shared taxa."""

    def __ilshift__(self, other):
        """c <<= d add counts of d but only in existing items in c."""
        if isinstance(self, col.Counter) and isinstance(other, col.Counter):
            for counter in other:
                if counter in self:
                    self[counter] += other[counter]
            return self
        return NotImplemented

    def __and__(self, other):
        """c & d add counts only for existing items in both c & d."""
        if isinstance(self, col.Counter) and isinstance(other, col.Counter):
            result: SharedCounter = SharedCounter()
            for item in other:
                if item in self:
                    result[item] = self[item] + other[item]
            return result
        return NotImplemented

    def __iand__(self, other):
        """c &= d add counts only for existing items in both c & d."""
        self = SharedCounter.__and__(self, other)
        return self

    def __floordiv__(self, other):
        """c // i floor divide each element of c by integer i."""
        if isinstance(self, col.Counter) and isinstance(other, int):
            result: SharedCounter = SharedCounter()
            for item in self:
                result[item] = self[item] // other
            return self
        return NotImplemented

    def __rfloordiv__(self, other):
        return SharedCounter.__floordiv__(self, other)

    def __ifloordiv__(self, other):
        """c //= i floor divide each element of c by integer i."""
        if isinstance(self, col.Counter) and isinstance(other, int):
            for counter in self:
                self[counter] //= other
            return self
        return NotImplemented


def read_report(report_file: str) -> Tuple[str, Counter[TaxId],
                                           Dict[TaxId, Rank]]:
    """
    Read Centrifuge/Kraken report file

    Args:
        report_file: report file name

    Returns:
        log string, abundances counter, taxlevel dict

    """
    output: io.StringIO = io.StringIO(newline='')
    abundances: Counter[TaxId] = col.Counter()
    level_dic = {}
    output.write(f'\033[90mLoading report file {report_file}...\033[0m')
    try:
        with open(report_file, 'r') as file:
            for report_line in file:
                _, _, taxnum, taxlev, _tid, _ = report_line.split('\t')
                tid = TaxId(_tid)
                abundances[tid] = int(taxnum)
                level_dic[tid] = Rank.centrifuge(taxlev)
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        report_file + '"')
    else:
        output.write('\033[92m OK! \033[0m\n')
    return output.getvalue(), abundances, level_dic


def write_lineage(parents: Parents,
                  names: Dict[TaxId, str],
                  tree: TaxTree,
                  lineage_file: str,
                  nodes: Counter[TaxId],
                  collapse: bool = True) -> str:
    """
    Writes a lineage file understandable by Krona.

    Args:
        parents: dictionary of parents for every TaxId.
        names: dictionary of names for every TaxId.
        tree: a TaxTree structure.
        lineage_file: name of the lineage file
        nodes: a counter for TaxIds
        collapse: This bool controls the collapse of taxid 131567
            (cellular organisms) and is True by default

    Returns: A string with the output messages

    """
    log, taxids_dic = tree.get_lineage(parents, iter(nodes))
    output: io.StringIO = io.StringIO(newline='')
    output.write(log)
    if collapse:  # Collapse taxid 131567 (cellular organisms) if desired
        for tid in taxids_dic:
            if len(taxids_dic[tid]) > 2:  # Not collapse for unclassified
                try:
                    taxids_dic[tid].remove(
                        TaxId('131567'))  # TaxId cellular organisms
                except ValueError:
                    pass
    lineage_dic = {tax_id: [names[tid] for tid in taxids_dic[tax_id]]
                   for tax_id in taxids_dic}
    output.write(f'  \033[90mSaving lineage file {lineage_file} with '
                 f'{len(nodes)} nodes...\033[0m')
    with open(lineage_file, 'w', newline='') as tsv_handle:
        tsvwriter = csv.writer(tsv_handle, dialect='krona')
        tsvwriter.writerow(["#taxID", "Lineage"])
        for tid in nodes:
            counts: str = str(nodes[tid])  # nodes[tid] is a int
            row: List[Union[TaxId, str]] = [counts, ]
            row.extend(lineage_dic[tid])
            tsvwriter.writerow(row)
    output.write('\033[92m OK! \033[0m\n')
    return output.getvalue()


def process_report(*args,
                   **kwargs
                   ) -> Tuple[Filename, TaxTree, TaxLevels, Counter[TaxId]]:
    """
    Process Centrifuge/Kraken report files (to be usually called in parallel!).
    """
    # Recover input and parameters
    filerep: Filename = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    parents: Parents = taxonomy.parents
    names: Names = taxonomy.names
    mintaxa: int = kwargs['mintaxa']
    collapse: bool = taxonomy.collapse
    including: Set[TaxId] = taxonomy.including
    excluding: Set[TaxId] = taxonomy.excluding
    verb: bool = kwargs['verb']
    output: io.StringIO = io.StringIO(newline='')

    # Read Centrifuge/Kraken report file to get abundances
    log, abundances, _ = read_report(filerep)
    output.write(log)

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy, abundances, _ROOT, [])  # Grow tax tree from root node
    output.write('\033[92m OK! \033[0m\n')

    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    tree.prune(mintaxa, None, collapse, verb)
    output.write('\033[92m OK! \033[0m\n')

    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    new_abund: Counter[TaxId] = col.Counter()
    ranks: Ranks = Ranks({})
    tree.get_taxa(new_abund, ranks,
                  mindepth=0, maxdepth=0,
                  include=including,  # ('2759',),  # ('135613',),
                  exclude=excluding)  # ('9606',))  # ('255526',))
    new_abund = +new_abund  # remove zero and negative counts
    taxlevels: TaxLevels = Rank.ranks_to_taxlevels(ranks)
    output.write('\033[92m OK! \033[0m\n')

    # Write the lineage file
    filename = Filename(filerep + _LINEAGEXT)
    log = write_lineage(parents, names, tree, filename, new_abund, collapse)
    output.write(log)
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return filename, tree, taxlevels, new_abund


def process_rank(*args, **kwargs):
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """
    # Recover input and parameters
    rank: Rank = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    parents = taxonomy.parents
    names = taxonomy.names
    collapse = taxonomy.collapse
    including = taxonomy.including
    excluding = taxonomy.excluding
    mintaxa = kwargs['mintaxa']
    trees: Dict[Sample, TaxTree] = kwargs['trees']
    taxids: Dict[Sample, TaxLevels] = kwargs['taxids']
    reports: List[Sample] = kwargs['reports']
    control: bool = kwargs['control']

    # Declare/define variables
    control_abundance: Counter[TaxId]
    #ranks: Ranks
    filenames: List[Filename] = []
    report: Sample
    filename: Filename
    log: str
    shared_abundance: SharedCounter = SharedCounter()
    shared_ctrl_abundance: SharedCounter = SharedCounter()
    output: io.StringIO = io.StringIO(newline='')

    # Cross analysis iterating by report: exclusive and part of shared&ctrl
    for i, report in zip(range(len(reports)), reports):
        exclude: Set[TaxId] = set()
        # Get taxids at this rank that are present in the other samples
        for sample in (_report for _report in reports if _report != report):
            exclude.update(taxids[sample][rank])
        exclude.update(excluding)  # Add explicit excluding taxa if any
        output.write(f'  \033[90mExcluding {len(exclude)} taxa from '
                     f'{report} at rank "{rank.name.lower()}"\n'
                     '  \033[90mFiltering shared taxa...\033[0m')
        leveled_tree = copy.deepcopy(trees[report])
        leveled_tree.prune(mintaxa, rank)  # Prune the copy to the desired rank
        # Get abundance for the exclusive analysis
        exclusive_abundance: Counter[TaxId] = col.Counter()
        leveled_tree.get_taxa(exclusive_abundance,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=exclude,  # Extended exclusion set
                              just_level=rank,
                              )
        # Get partial abundance for the shared analysis
        sub_shared_abundance: SharedCounter = SharedCounter()
        leveled_tree.get_taxa(sub_shared_abundance,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=excluding,
                              just_level=rank,
                              )
        output.write('\033[92m OK! \033[0m\n')
        exclusive_abundance = +exclusive_abundance  # remove counts <= 0
        filename = Filename(f'EXCLUSIVE_{rank.name.lower()}_'
                            f'{report}{_LINEAGEXT}')
        log = write_lineage(parents, names, trees[report],
                            filename, exclusive_abundance, collapse)
        output.write(log)
        filenames.append(filename)
        # Shared and shared-control taxa partial evaluations
        if i == 0:  # 1st iteration: Initialize shared_abundance
            shared_abundance.update(sub_shared_abundance)
        elif i == 1:  # 2nd iteration: Initialize shared_ctrl_abundance
            shared_abundance &= sub_shared_abundance
            shared_ctrl_abundance.update(sub_shared_abundance)
        else:  # Accumulate shared abundances
            shared_abundance &= sub_shared_abundance
            shared_ctrl_abundance &= sub_shared_abundance

    # Shared taxa final analysis
    output.write('  \033[90mAnalyzing shared taxa...\033[0m')
    shared_abundance //= len(reports)  # Get averaged abundance
    output.write('\033[92m OK! \033[0m\n')
    if shared_abundance:
        output.write(f'  \033[90mIncluding {len(shared_abundance)} taxa at '
                     f'rank "{rank.name.lower()}"\033[0m\n'
                     '  \033[90mBuilding taxonomy tree...\033[0m')
        tree = TaxTree()
        tree.grow(taxonomy, shared_abundance)
        output.write('\033[92m OK! \033[0m\n')
        filename = Filename(f'SHARED_{rank.name.lower()}{_LINEAGEXT}')
        log = write_lineage(parents, names, tree,
                            filename, shared_abundance, collapse)
        output.write(log)
        filenames.append(filename)

    if control:
        # Control sample substraction
        output.write(f'  \033[90m{reports[0]} is selected as control '
                     f'sample for subtractions.\033[0m\n')
        # Get taxids at this rank that are present in the ctrl sample
        exclude = set(taxids[reports[0]][rank])
        exclude.update(excluding)  # Add explicit excluding taxa if any
        for report in reports[1::]:
            output.write(f'  \033[90mExcluding {len(exclude)} ctrl taxa '
                         f'from {report} at rank "{rank.name.lower()}"\n'
                         '  \033[90mFiltering taxa from control...\033[0m')
            leveled_tree = copy.deepcopy(trees[report])
            leveled_tree.prune(mintaxa, rank)
            control_abundance = col.Counter()
            leveled_tree.get_taxa(control_abundance,
                                  include=including,
                                  exclude=exclude,  # Extended exclusion set
                                  just_level=rank,
                                  )
            output.write('\033[92m OK! \033[0m\n')
            control_abundance = +control_abundance  # remove counts <= 0
            filename = Filename(f'CONTROL_{rank.name.lower()}_'
                                f'{report}{_LINEAGEXT}')
            log = write_lineage(parents, names, trees[report],
                                filename, control_abundance, collapse)
            output.write(log)
            filenames.append(filename)

        # Shared-control taxa final analysis
        output.write('  \033[90mAnalyzing shared-control taxa...\033[0m')
        shared_ctrl_abundance //= (len(reports) - 1)  # Get averaged abundance
        output.write('\033[92m OK! \033[0m\n')
        if shared_ctrl_abundance:
            output.write(
                f'  \033[90mIncluding {len(shared_ctrl_abundance)} taxa at '
                f'rank "{rank.name.lower()}"\033[0m\n'
                '  \033[90mBuilding taxonomy tree...\033[0m')
            tree = TaxTree()
            tree.grow(taxonomy, shared_ctrl_abundance)
            tree.prune(mintaxa, rank)
            new_shared_ctrl_abundance: Counter[TaxId] = col.Counter()
            tree.get_taxa(new_shared_ctrl_abundance,
                          include=including,
                          exclude=exclude,  # Extended exclusion set
                          just_level=rank,
                          )
            output.write('\033[92m OK! \033[0m\n')
            filename = Filename(f'SHARED_CONTROL_{rank.name.lower()}'
                                f'{_LINEAGEXT}')
            log = write_lineage(parents, names, tree,
                                filename, new_shared_ctrl_abundance, collapse)
            output.write(log)
            filenames.append(filename)

    # Print output and return
    print(output.getvalue())
    sys.stdout.flush()
    return filenames


def write_krona(samples: List[Sample],
                outputs: Dict[Rank, Filename],
                htmlfile: Filename = Filename('Output' + _HTML_SUFFIX),
                ):
    """Generate the Krona html file calling ktImportText."""
    subprc = ["ktImportText"]
    subprc.extend(samples)
    try:
        subprc.extend([outputs[level][i]
                       for level in Rank.selected_ranks
                       for i in range(len(outputs[level]))])
    except KeyError:
        pass
    subprc.extend(["-o", htmlfile])
    try:
        subprocess.run(subprc, check=True)
    except subprocess.CalledProcessError:
        print('\n\033[91mERROR!\033[0m ktImportText: ' +
              'returned a non-zero exit status (Krona plot built failed)')


def main():
    """Main entry point to recentrifuge."""
    # Argument Parser Configuration
    parser = argparse.ArgumentParser(
        description='Post-process Centrifuge/Kraken output',
        epilog=f'%(prog)s  - {__author__} - {__date__}'
    )
    parser.add_argument(
        '-V', '--version',
        action='version',
        version=f'%(prog)s release {__version__} ({__date__})'
    )
    parser.add_argument(
        '-f', '--filerep',
        action='append',
        metavar='FILE',
        required=True,
        help=('Centrifuge/Kraken report files ' +
              '(multiple -f is available to include several samples in plot')
    )
    parser.add_argument(
        '-v', '--verbose',
        action='count',
        default=0,
        help='increase output verbosity'
    )
    parser.add_argument(
        '-n', '--nodespath',
        action='store',
        metavar='PATH',
        default='./',
        help=('path for the nodes information files (nodes.dmp and names.dmp' +
              ' from NCBI')
    )
    parser.add_argument(
        '-m', '--mintaxa',
        action='store',
        metavar='INT',
        default=_DEFMINTAXA,
        help=('minimum taxa to avoid collapsing one level to the parent one ' +
              ('(%i per default)' % _DEFMINTAXA))
    )
    parser.add_argument(
        '-k', '--nokollapse',
        action='store_true',
        help='Show the "cellular organisms" taxon (collapsed by default)'
    )
    parser.add_argument(
        '-o', '--outhtml',
        action='store',
        metavar='FILE',
        default=None,
        help='HTML Krona file (default: name inferred from report files)'
    )
    parser.add_argument(
        '-i', '--include',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to include a taxon and all underneath ' +
              '(multiple -i is available to include several taxid). ' +
              'By default all the taxa is considered for inclusion.')
    )
    parser.add_argument(
        '-x', '--exclude',
        action='append',
        metavar='TAXID',
        type=TaxId,
        default=[],
        help=('NCBI taxid code to exclude a taxon and all underneath ' +
              '(multiple -x is available to exclude several taxid)')
    )
    parser.add_argument(
        '-s', '--sequential',
        action='store_true',
        help='Deactivate parallel processing.'
    )
    parser.add_argument(
        '-a', '--avoidcross',
        action='store_true',
        help='Avoid cross analysis.'
    )
    parser.add_argument(
        '-c', '--control',
        action='store_true',
        help='Take the first sample as negative control.'
    )

    # Parse arguments
    args = parser.parse_args()
    reports = args.filerep
    verb = args.verbose
    nodesfile = os.path.join(args.nodespath, _NODESFILE)
    namesfile = os.path.join(args.nodespath, _NAMESFILE)
    mintaxa = int(args.mintaxa)
    collapse = not args.nokollapse
    excluding: Set[TaxId] = set(args.exclude)
    including: Set[TaxId] = set(args.include)
    sequential = args.sequential
    avoidcross = args.avoidcross
    control = args.control
    htmlfile: Filename = args.outhtml
    if not htmlfile:
        htmlfile = reports[0].split('_mhl')[0] + _HTML_SUFFIX

    # Program header and chdir
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Load NCBI nodes, names and build children
    ncbi: Taxonomy = Taxonomy(nodesfile, namesfile,
                              collapse, excluding, including)

    # Declare variables that will hold results for the samples analyzed
    trees: Dict[Sample, TaxTree] = {}
    abundances: Dict[Sample, Counter[TaxId]] = {}
    taxids: Dict[Sample, TaxLevels] = {}
    samples: List[Sample] = []
    # Processing of report files in parallel
    print('\033[90mPlease, wait, processing files in parallel...\033[0m\n')
    kwargs = {'taxonomy': ncbi, 'mintaxa': mintaxa, 'verb': verb}
    # Enable parallelization with 'spawn' under known platforms
    if platform.system() and not sequential:  # Only for known platforms
        mpctx = mp.get_context('spawn')  # Important for OSX&Win
        with mpctx.Pool(processes=min(os.cpu_count(),
                                      len(reports))) as pool:
            async_results = [pool.apply_async(
                process_report,
                args=[filerep],
                kwds=kwargs
            ) for filerep in reports]
            for sample, (filename, trees[sample], taxids[sample],
                         abundances[sample]) in zip(reports,
                                                    [r.get() for r in
                                                     async_results]):
                samples.append(filename)
    else:  # sequential processing of each sample
        for sample in reports:
            (filename,
             trees[sample],
             taxids[sample],
             abundances[sample]) = process_report(sample, **kwargs)
            samples.append(filename)

    # Cross analysis of samples in parallel by taxlevel
    outputs: Dict[Rank, Filename] = {}  # dict: rank to output files (w/o ext)
    # Avoid if just a single report file of explicitly stated by flag
    if len(reports) > 1 and not avoidcross:
        print('\033[90mPlease, wait. ' +
              'Performing cross analysis in parallel...\033[0m\n')
        kwargs.update({'trees': trees, 'taxids': taxids,
                       'abundances': abundances, 'reports': reports,
                       'control': control})
        if platform.system() and not sequential:  # Only for known platforms
            mpctx = mp.get_context('spawn')  # Important for OSX&Win
            with mpctx.Pool(processes=min(os.cpu_count(), len(
                    Rank.selected_ranks))) as pool:
                async_results = [pool.apply_async(
                    process_rank,
                    args=[level],
                    kwds=kwargs
                ) for level in Rank.selected_ranks]
                for level, (outputs[level]) in zip(
                        Rank.selected_ranks,
                        [r.get() for r in async_results]):
                    pass
        else:  # sequential processing of each selected rank
            for level in Rank.selected_ranks:
                outputs[level] = process_rank(level, **kwargs)

    # Generate Krona plot with all the results
    write_krona(samples, outputs, htmlfile)


if __name__ == '__main__':
    main()
