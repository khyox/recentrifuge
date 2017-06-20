#!/usr/bin/env python3
"""
Post-process Centrifuge/Kraken report files for Krona.
"""
# pylint: disable=not-an-iterable, bad-reversed-sequence, E0203
__version__ = '0.6.2'
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
from typing import NewType, Iterator, Union

# Type annotations
# pylint: disable=invalid-name
TaxId = NewType('TaxId', str)
Parents = NewType('Parents', Dict[TaxId, TaxId])
Ranks = NewType('Ranks', Dict[TaxId, 'Rank'])
Names = NewType('Names', Dict[TaxId, str])
Children = NewType('Children', Dict[TaxId, Dict[TaxId, int]])
# pylint: enable=invalid-name

# Predefined internal constants
_PATH = '.'
_NODESFILE = 'nodes.dmp'
_NAMESFILE = 'names.dmp'
_LINEAGEXT = '.lineage'
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
    def invert_dict(cls, ranks: dict) -> dict:
        """Invert dictionary from items->Rank to Rank->set(items)"""
        return {cls[taxlevel]: {item for item in ranks if
                                ranks[item] is cls[taxlevel]} for
                taxlevel in cls.__members__}  # type: ignore

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
                 excluding: Tuple[TaxId, ...],
                 including: Tuple[TaxId, ...],
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
            including = (_ROOT,)
        self.including: Tuple[TaxId, ...] = including
        if excluding:
            print('List of taxa (and below) to be excluded:')
            print('\t\tTaxId\tScientific Name')
            for taxid in excluding:
                print(f'\t\t{taxid}\t{self.names[taxid]}')
        self.excluding: Tuple[TaxId, ...] = excluding

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
             abundances: Counter[TaxId],
             taxid: TaxId,
             path: List[TaxId]
             ) -> None:
        """Recursive function to build a taxonomy tree"""
        if taxid not in path:  # Avoid loops for repeated taxid (like root)
            self[taxid] = TaxTree(counts=abundances.get(taxid, 0),
                                  rank=taxonomy.get_rank(taxid))
            if taxid in taxonomy.children:  # taxid has children
                for child in taxonomy.children[taxid]:
                    self[taxid].grow(taxonomy, abundances, child, path+[taxid])

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
                 abundance: Counter,
                 taxlevels: Ranks,
                 mindepth: int = 0,
                 maxdepth: int = 0,
                 include: Tuple = (),
                 exclude: Tuple = (),
                 just_level: Rank = None,
                 _in_branch: bool = False
                 ) -> None:
        """
        Recursively get the taxa between min and max depth levels.

        Args:
            abundance: Input/output Counter; at 1st entry it's empty.
            taxlevels: Input/output dict; at 1st entry it's empty.
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
                    taxlevels[tid] = self[tid].taxlevel
                if self[tid]:
                    self[tid].get_taxa(abundance, taxlevels,
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


def process_report(*args, **kwargs):
    """
    Process Centrifuge/Kraken report files (to be usually called in parallel!).
    """
    # Recover input and parameters
    filerep = args[0]
    taxonomy = kwargs['taxonomy']
    parents = taxonomy.parents
    names = taxonomy.names
    mintaxa = kwargs['mintaxa']
    collapse = taxonomy.collapse
    including = taxonomy.including
    excluding = taxonomy.excluding
    verb = kwargs['verb']
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
    taxlevels: Ranks = Ranks({})
    tree.get_taxa(new_abund, taxlevels,
             mindepth=0, maxdepth=0,
             include=including,  # ('2759',),  # ('135613',),
             exclude=excluding)  # ('9606',))  # ('255526',))
    new_abund = +new_abund  # remove zero and negative counts
    taxid: Dict[Rank, TaxId] = Rank.invert_dict(taxlevels)
    output.write('\033[92m OK! \033[0m\n')

    # Write the lineage file
    filename = filerep + _LINEAGEXT
    log = write_lineage(parents, names, tree, filename, new_abund, collapse)
    output.write(log)
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return filename, tree, taxid, new_abund


def process_rank(*args, **kwargs):
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """
    # Recover input and parameters
    rank: Rank = args[0]
    taxonomy = kwargs['taxonomy']
    parents = taxonomy.parents
    names = taxonomy.names
    collapse = taxonomy.collapse
    including = taxonomy.including
    excluding = taxonomy.excluding
    mintaxa = kwargs['mintaxa']
    tree_dic = kwargs['tree_dic']
    taxid_dic = kwargs['taxid_dic']
    abund_dic = kwargs['abund_dic']
    filerep_lst = kwargs['filerep_lst']
    verb = kwargs['verb']
    abundances: Counter[TaxId]
    ranks: Ranks

    filen_lst = []
    filename: str
    output: io.StringIO = io.StringIO(newline='')
    # Exclusive (not shared) taxa analysis
    for filerep in filerep_lst:
        exclude_set: set = set()
        for filename in (fn for fn in filerep_lst if fn != filerep):
            for subrank in rank.ranks_from_specific:
                exclude_set.update(taxid_dic[filename][subrank])
        exclude_set.update(excluding)
        exclude_tup: Tuple[TaxId, ...] = tuple(exclude_set)
        if exclude_tup:
            output.write(f'  \033[90mExcluding {len(exclude_tup)} taxa from '
                         f'{filerep} at rank "{rank.name.lower()}"\n'
                         '  \033[90mFiltering shared taxa...\033[0m')
        abundances = col.Counter()
        ranks = Ranks({})
        leveled_tree = copy.deepcopy(tree_dic[filerep])
        leveled_tree.prune(mintaxa, rank)
        leveled_tree.get_taxa(abundances, ranks,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=exclude_tup,
                              just_level=rank,
                 )
        output.write('\033[92m OK! \033[0m\n')
        abundances = +abundances  # remove zero and negative counts
        filename = f'EXCLUSIVE_{rank.name.lower()}_{filerep}{_LINEAGEXT}'
        log: str = write_lineage(parents, names, tree_dic[filerep],
                                 filename, abundances, collapse)
        output.write(log)
        filen_lst.append(filename)

    # Shared (in common) taxa analysis
    output.write('  \033[90mAnalyzing shared taxa...\033[0m')
    include: Set[TaxId] = set()
    for subrank in rank.ranks_from_general:
        include_subset = taxid_dic[filerep_lst[0]][subrank]
        for filerep in filerep_lst:
            include_subset &= taxid_dic[filerep][subrank]  # Acc intersection
        include |= include_subset

    # weight = 1.0/len(filerep_lst)
    # Calculate shared taxa abundance as mean of the abs counts in each dataset
    abundances = col.Counter({taxid: sum([abund_dic[fn][taxid]
                                        for fn in filerep_lst]) // len(
        filerep_lst) for taxid in include})
    #    for filerep in filerep_lst:
    #        for tid in include_set:
    #            abun_cnt[tid] += weight*abund_dic[filerep][tid]
    output.write('\033[92m OK! \033[0m\n')
    if abundances:
        output.write(f'  \033[90mIncluding {len(abundances)} taxa at '
                     f'rank "{rank.name.lower()}"\033[0m\n'
                     '  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(taxonomy, abundances, _ROOT, [])
    output.write('\033[92m OK! \033[0m\n')
    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    tree.prune(mintaxa, None, collapse, verb)
    output.write('\033[92m OK! \033[0m\n')
    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    abundances = col.Counter()
    ranks = Ranks({})
    tree.get_taxa(abundances, ranks,
                  mindepth=0, maxdepth=0,
                  include=including, exclude=excluding,
                  just_level=rank,
                  )
    abundances = +abundances  # remove zero and negative counts
    output.write('\033[92m OK! \033[0m\n')
    # Save the taxonomy in tsv file
    filename = f'SHARED_{rank.name.lower()}{_LINEAGEXT}'
    log = write_lineage(parents, names, tree,
                        filename, abundances, collapse)
    output.write(log)
    filen_lst.append(filename)
    print(output.getvalue())
    sys.stdout.flush()

    # Return
    return filen_lst


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

    # Parse arguments
    args = parser.parse_args()
    filerep_lst = args.filerep
    verb = args.verbose
    nodesfile = os.path.join(args.nodespath, _NODESFILE)
    namesfile = os.path.join(args.nodespath, _NAMESFILE)
    mintaxa = int(args.mintaxa)
    collapse = not args.nokollapse
    excluding: Tuple[TaxId, ...] = tuple(args.exclude)
    including: Tuple[TaxId, ...] = tuple(args.include)
    sequential = args.sequential
    avoidcross = args.avoidcross
    htmlfile = args.outhtml
    if not htmlfile:
        htmlfile = filerep_lst[0].split('_mhl')[0] + '.cf2krona.html'

    # Program header and chdir
    print(f'\n=-= {sys.argv[0]} =-= v{__version__} =-= {__date__} =-=\n')
    sys.stdout.flush()

    # Load NCBI nodes, names and build children
    ncbi: Taxonomy = Taxonomy(nodesfile, namesfile,
                              collapse, excluding, including)

    tree_dic = {}  # dict with the tree of every sample
    abund_dic = {}  # sample dict with col.Counter with abundances of taxa
    #    level_dic = {}  # sample dict with the taxonomic levels of taxa
    taxid_dic = {}  # sample dict with dict of taxids for every taxlevel
    filen_lst = []  # sample dict with the name of the output file (w/o ext)
    # Processing of report files in parallel
    print('\033[90mPlease, wait, processing files in parallel...\033[0m\n')
    kwargs = {'taxonomy': ncbi, 'mintaxa': mintaxa, 'verb': verb}
    # Enable parallelization with 'spawn' under known platforms
    if platform.system() and not sequential:  # Only for known platforms
        mpctx = mp.get_context('spawn')  # Important for OSX&Win
        with mpctx.Pool(processes=min(os.cpu_count(),
                                      len(filerep_lst))) as pool:
            async_results = [pool.apply_async(
                process_report,
                args=[filerep],
                kwds=kwargs
            ) for filerep in filerep_lst]
            for filerep, (filen, tree_dic[filerep], taxid_dic[filerep],
                          abund_dic[filerep]) in zip(filerep_lst,
                                                     [r.get() for r in
                                                      async_results]):
                filen_lst.append(filen)
    else:
        for filerep in filerep_lst:
            filen, tree_dic[filerep], taxid_dic[filerep], abund_dic[filerep] = \
                process_report(filerep, **kwargs)
            filen_lst.append(filen)

    # Cross analysis of samples in parallel by taxlevel
    filen_diff_dic = {}  # taxlevel dict with the name of output file (w/o ext)
    # Avoid if just a single report file of explicitly stated by flag
    if len(filerep_lst) > 1 and not avoidcross:
        print('\033[90mPlease, wait. ' +
              'Performing cross analysis in parallel...\033[0m\n')
        kwargs.update({'tree_dic': tree_dic, 'taxid_dic': taxid_dic,
                       'abund_dic': abund_dic, 'filerep_lst': filerep_lst})
        if platform.system() and not sequential:  # Only for known platforms
            mpctx = mp.get_context('spawn')  # Important for OSX&Win
            with mpctx.Pool(processes=min(os.cpu_count(), len(
                    Rank.selected_ranks))) as pool:
                async_results = [pool.apply_async(
                    process_rank,
                    args=[level],
                    kwds=kwargs
                ) for level in Rank.selected_ranks]
                for level, (filen_diff_dic[level]) in zip(
                        Rank.selected_ranks,
                        [r.get() for r in async_results]):
                    pass
        else:
            for level in Rank.selected_ranks:
                filen_diff_dic[level] = process_rank(level, **kwargs)

    # Generate the Krona html file calling ktImportText
    subprc = ["ktImportText"]
    subprc.extend(filen_lst)
    try:
        subprc.extend([filen_diff_dic[level][i]
                       for level in Rank.selected_ranks
                       for i in range(len(filen_diff_dic[level]))])
    except KeyError:
        pass
    subprc.extend(["-o", htmlfile])
    try:
        subprocess.run(subprc, check=True)
    except subprocess.CalledProcessError:
        print('\n\033[91mERROR!\033[0m ktImportText: ' +
              'returned a non-zero exit status (Krona plot built failed)')


if __name__ == '__main__':
    # for name, member in Rank.__members__.items():
    #     print(name, member, member.value)
    main()
