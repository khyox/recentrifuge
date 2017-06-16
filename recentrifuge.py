#!/usr/bin/env python3
"""
Post-process Centrifuge/Kraken report files for Krona.
"""
# pylint: disable=not-an-iterable, bad-reversed-sequence

__version__ = '0.6.0'
__author__ = 'Jose Manuel Marti'
__date__ = 'Jun 2017'

import argparse
import collections as col
import csv
import io
import multiprocessing as mp
import os
import platform
import subprocess
import sys
from enum import Enum, unique
from typing import Iterable, Dict, Counter, Tuple, List, Any
from typing import NewType, Iterator, Union

# Type annotations
TaxId = NewType('TaxId', str)  # pylint: disable=invalid-name

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


@unique
class TaxLevel(Enum):
    """Enumeration with codes for taxonomical levels."""
    # pylint: disable=invalid-name
    U = 'Unclassified'
    D = 'Domain'
    K = 'Kingdom'
    P = 'Phylum'
    C = 'Class'
    O = 'Order'
    F = 'Family'
    G = 'Genus'
    S = 'Species'
    H = 'otHer'
    W = 'unknoWn'

    # pylint: enable=invalid-name


    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)

    @classmethod
    def centrifuge(cls, tax_level: str) -> 'TaxLevel':
        """Transforms Centrifuge codes for taxonomical levels"""
        taxonomic_level: TaxLevel
        if tax_level == '-':
            taxonomic_level = cls.H
        else:
            try:
                taxonomic_level = cls[tax_level]
            except KeyError:
                raise UnsupportedTaxLevelError(
                    f'Unknown tax level {tax_level}')
        return taxonomic_level

    @classmethod
    def invert_dict(cls, tax_levels: dict) -> dict:
        """Invert dictionary from items->TaxLevel to TaxLevel->set(items)"""
        return {cls[taxlevel]: {item for item in tax_levels if
                                tax_levels[item] is cls[taxlevel]} for
                taxlevel in cls.__members__}  # type: ignore

    @classproperty
    def selected_taxlevels(cls):  # pylint: disable=no-self-argument
        """Tax levels selected for deep analysis and comparisons"""
        return [cls.S, cls.G, cls.F, cls.O, cls.C, cls.P, cls.K, cls.D]

    @property
    def taxlevels_from_specific(self) -> Iterator['TaxLevel']:
        """Generator returning selected taxlevels from specific to general."""
        for taxlevel in self.__class__.selected_taxlevels:  # type: ignore
            yield taxlevel
            if taxlevel is self:
                break

    @property
    def taxlevels_from_general(self) -> Iterator['TaxLevel']:
        """Generator returning selected taxlevels from general to specific."""
        for taxlevel in reversed(
                self.__class__.selected_taxlevels):  # type: ignore
            yield taxlevel
            if taxlevel is self:
                break


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
        self.parents: Dict[TaxId, TaxId] = {}
        self.names: Dict[TaxId, str] = {}
        self.children: Dict[TaxId, Dict[TaxId, int]] = {}
        self.collapse: bool = collapse

        # Initialization methods
        self.read_parents(nodes_file)
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

    def read_parents(self, nodes_file: str) -> None:
        """Build dict with parent for a given taxid (key)"""
        print('\033[90mLoading NCBI nodes...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(nodes_file, 'r') as file:
                for line in file:
                    tid, parent, *_ = line.split('\t|\t')
                    self.parents[tid] = parent
        except:
            raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                            nodes_file + '"')
        else:
            print('\033[92m OK! \033[0m')

    def read_names(self, names_file: str) -> None:
        """Build dict with name for a given taxid (key)"""
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
        """Build dict of children for a given parent taxid (key)"""
        print('\033[90mBuilding dict of parent to children taxa...\033[0m',
              end='')
        sys.stdout.flush()
        self.children: Dict[TaxId, Dict[TaxId, int]] = {}
        for tid in self.parents:
            if self.parents[tid] not in self.children:
                self.children[self.parents[tid]] = {}
            self.children[self.parents[tid]][tid] = 0
        print('\033[92m OK! \033[0m')


class TaxTree(dict):
    """Nodes of a taxonomical tree"""

    def __init__(self, *args,
                 counts: int = 0,
                 taxlevel: TaxLevel = TaxLevel.W
                 ) -> None:
        super().__init__(args)
        self.counts: int = counts
        self.taxlevel: TaxLevel = taxlevel

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

    def grow(self, children, parentid, path, abundances, taxlevels):
        """Recursive function to build the taxonomy tree"""
        if parentid not in path:  # Avoid loops for repeated taxid (like root)
            self[parentid] = TaxTree(counts=abundances.get(parentid, 0),
                                     taxlevel=taxlevels.get(parentid,
                                                            TaxLevel.W))
            if parentid in children:
                for child in children[parentid]:
                    self[parentid].grow(children, child, path + [parentid],
                                        abundances, taxlevels)


def prune_tree(subtree, nmin=1, collapse=False, verb=False):
    """Recursive function to collapse low abundant taxa of the taxonomy tree"""
    for tid in list(subtree):
        if subtree[tid]:  # Check if every subtree tid is already a leaf
            if prune_tree(subtree[tid],
                          nmin=nmin) and subtree[tid].counts < nmin:
                if collapse:
                    subtree.counts += subtree[
                        tid].counts  # Acc abundance in higher tax
                subtree.pop(tid)  # Prune branch if it has no leafs anymore
            elif verb:
                print("NOT pruning branch", tid, "with", subtree[tid].counts)
        else:
            if subtree[tid].counts < nmin:
                if collapse:
                    subtree.counts += subtree[
                        tid].counts  # Accumulate abundance in higher tax
                subtree.pop(tid)  # Prune leaf
            elif verb:
                print("NOT pruning leaf", tid, "with", subtree[tid].counts)
    return not subtree  # returns True if this subtree is a leaf


def get_taxa(subtree: TaxTree,
             abundance: Counter,
             taxlevels: Dict[TaxId, TaxLevel],
             mindepth: int = 0,
             maxdepth: int = 0,
             include: Tuple = (),
             exclude: Tuple = (),
             just_level: TaxLevel = None,
             _in_branch: bool = False
             ) -> None:
    """
    Recursive function to get the taxa between min and max depth levels

    Args:
        subtree: TaxTree structure.
        abundance: Is a input/output Counter: at entry it's empty.
        taxlevels: Is a input/output dict: at entry it's empty.
        mindepth: 0 gets the taxa from the very beginning depth level.
        maxdepth: 0 does not stop the search at any depth level.
        include: contains the root taxid of the subtrees to be
            included. If it is empty (default) all the taxa is included
            (except explicitly excluded).
        exclude: contains the root taxid of the subtrees to be excluded
        just_level: If set, just taxa in this taxlevel will be counted.
        _in_branch: is like a static variable to tell the recursive
            function that it is in a subtree that is a branch of a
            taxon in the include list.

    Returns: None

    """
    mindepth -= 1
    if maxdepth != 1:
        maxdepth -= 1
        for tid in subtree:
            in_branch: bool = ((_in_branch or  # called with flag activated? or
                                not include or  # include by default? or
                                (tid in include))  # tid is to be included?
                               and tid not in exclude)  # and not to exclude
            if (mindepth <= 0
                and in_branch
                and (just_level is None
                     or subtree[tid].taxlevel is just_level)):
                abundance[tid] = subtree[tid].counts
                taxlevels[tid] = subtree[tid].taxlevel
            if subtree[tid]:
                get_taxa(subtree[tid], abundance, taxlevels,
                         mindepth, maxdepth,
                         include, exclude,
                         just_level, in_branch)


def read_report(report_file: str) -> Tuple[str, Counter,
                                           Dict[TaxId, TaxLevel]]:
    """
    Read Centrifuge/Kraken report file

    Args:
        report_file: report file name

    Returns: log string, abundances counter, taxlevel dict

    """
    output: io.StringIO = io.StringIO(newline='')
    abundances: Counter[TaxId] = col.Counter()
    level_dic = {}
    output.write(f'\033[90mLoading report file {report_file}...\033[0m')
    try:
        with open(report_file, 'r') as file:
            for report_line in file:
                _, _, taxnum, taxlev, tid, *_ = report_line.split('\t')
                abundances[tid] = int(taxnum)
                level_dic[tid] = TaxLevel.centrifuge(taxlev)
    except:
        raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                        report_file + '"')
    else:
        output.write('\033[92m OK! \033[0m\n')
    return output.getvalue(), abundances, level_dic


def trace_node(tree: TaxTree,
               targets: List[TaxId],
               nodes: List[TaxId]):
    """Recursive function to build a list of nodes from root to target taxid"""
    for tid in tree:
        if tree[tid]:
            nodes.append(tid)
            if targets[0] in tree[tid]:
                target = targets.pop()
                if target != _ROOT:  # Avoid to append twice the root node
                    nodes.append(target)
                break
            trace_node(tree[tid], targets, nodes)
            if not targets:
                break
            else:
                nodes.pop()


def get_lineage(parents: Dict[TaxId, TaxId],
                tree: TaxTree,
                tid_iter: Iterable):
    """
    Build dict with taxid as keys, whose values are the list of nodes
        in the tree in the path from root to such a taxid
    """
    output = io.StringIO('  \033[90mGetting lineage of taxa...\033[0m',
                         newline='')
    tid_node_dic: Dict[TaxId, List[TaxId]] = {}
    for tid in tid_iter:
        if tid == _ROOT:
            tid_node_dic[tid] = [_ROOT]
        elif tid in parents:
            node_lst: List[TaxId] = []
            target_lst = [tid]
            trace_node(tree, target_lst, node_lst)
            if node_lst:
                tid_node_dic[tid] = node_lst
            else:
                output.write('[\033[93mWARNING\033[0m: Failed trace_node '
                             f'of tid {tid}]\n')
        else:
            output.write('[\033[93mWARNING\033[0m: Discarded unknown tid '
                         f'{tid}]\n')
    output.write('\033[92m OK! \033[0m\n')
    return output.getvalue(), tid_node_dic


def write_lineage(parents: Dict[TaxId, TaxId],
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
    log, taxids_dic = get_lineage(parents, tree, iter(nodes))
    output = io.StringIO(log, newline='')
    if collapse:  # Collapse taxid 131567 (cellular organisms) if desired
        for tid in taxids_dic:
            if len(taxids_dic[tid]) > 2:  # Not collapse for unclassified
                try:
                    taxids_dic[tid].remove(
                        '131567')  # TaxId cellular organisms
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
    parents = kwargs['taxonomy'].parents
    children = kwargs['taxonomy'].children
    names = kwargs['taxonomy'].names
    mintaxa = kwargs['mintaxa']
    collapse = kwargs['taxonomy'].collapse
    including = kwargs['taxonomy'].including
    excluding = kwargs['taxonomy'].excluding
    verb = kwargs['verb']
    output: io.StringIO = io.StringIO(newline='')

    # Read Centrifuge/Kraken report file to get abundances
    log, abund, level = read_report(filerep)
    output.write(log)

    # Build taxonomy tree
    output.write('  \033[90mBuilding taxonomy tree...\033[0m')
    # tree_dic[filerep] = TaxTree([_ROOT, TaxTree()])
    tree = TaxTree()
    # grow_tree(_ROOT, [], abund_dic[filerep], tree_dic[filerep][_ROOT])
    tree.grow(children, _ROOT, [], abund, level)
    output.write('\033[92m OK! \033[0m\n')
    # print_tree(tree_dic[filerep])

    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    prune_tree(tree, nmin=mintaxa, collapse=collapse, verb=verb)
    output.write('\033[92m OK! \033[0m\n')

    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    abundances: Counter[TaxId] = col.Counter()
    taxlevels: Dict[TaxId, TaxLevel] = {}
    get_taxa(tree, abundances, taxlevels,
             mindepth=0, maxdepth=0,
             include=including,  # ('2759',),  # ('135613',),
             exclude=excluding)  # ('9606',))  # ('255526',))
    abundances = +abundances  # remove zero and negative counts
    taxid: Dict[TaxLevel, TaxId] = TaxLevel.invert_dict(taxlevels)
    output.write('\033[92m OK! \033[0m\n')

    # Write the lineage file
    filename = filerep + _LINEAGEXT
    log = write_lineage(parents, names, tree, filename, abundances, collapse)
    output.write(log)
    print(output.getvalue())
    sys.stdout.flush()
    # Return
    return filename, tree, taxid, abundances


def process_taxlevel(*args, **kwargs):
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """
    # Recover input and parameters
    level: TaxLevel = args[0]
    parents = kwargs['taxonomy'].parents
    children = kwargs['taxonomy'].children
    names = kwargs['taxonomy'].names
    mintaxa = kwargs['mintaxa']
    kollapse = kwargs['taxonomy'].collapse
    including = kwargs['taxonomy'].including
    excluding = kwargs['taxonomy'].excluding
    tree_dic = kwargs['tree_dic']
    taxid_dic = kwargs['taxid_dic']
    abund_dic = kwargs['abund_dic']
    filerep_lst = kwargs['filerep_lst']

    filen_lst = []
    output: io.StringIO = io.StringIO(newline='')
    # Exclusive (not shared) taxa analysis
    for filerep in filerep_lst:
        exclude_set: set = set()
        for filename in (fn for fn in filerep_lst if fn != filerep):
            for sublevel in level.taxlevels_from_specific:
                exclude_set.update(taxid_dic[filename][sublevel])
        exclude_set.update(excluding)
        exclude_tup: Tuple[TaxId, ...] = tuple(exclude_set)
        if exclude_tup:
            output.write(f'  \033[90mExcluding {len(exclude_tup)} taxa from '
                         f'{filerep} at taxlevel "{level.value}"\n'
                         '  \033[90mFiltering shared taxa...\033[0m')
        abundances: Counter[TaxId] = col.Counter()
        taxlevels: Dict[TaxId, TaxLevel] = {}
        get_taxa(tree_dic[filerep], abundances, taxlevels,
                 mindepth=0, maxdepth=0,
                 include=including,
                 exclude=exclude_tup,
                 just_level=level,
                 )
        output.write('\033[92m OK! \033[0m\n')
        abundances = +abundances  # remove zero and negative counts
        filen_diff = ('EXCLUSIVE_' + level.value + '_' + filerep + _LINEAGEXT)
        log: str = write_lineage(parents, names, tree_dic[filerep],
                                 filen_diff, abundances, kollapse)
        output.write(log)
        filen_lst.append(filen_diff)

    # Shared (in common) taxa analysis
    output.write('  \033[90mAnalyzing shared taxa...\033[0m')
    include_set = set()
    for sublevel in level.taxlevels_from_general:
        include_subset = taxid_dic[filerep_lst[0]][sublevel]
        for filerep in filerep_lst[1:]:
            include_subset &= taxid_dic[filerep][sublevel]
        include_set |= include_subset

    # weight = 1.0/len(filerep_lst)
    # Calculate shared taxa abundance as mean of the abs counts in each dataset
    abundances = col.Counter({taxid: sum([abund_dic[fn][taxid]
                                        for fn in filerep_lst]) // len(
        filerep_lst) for taxid in include_set})
    #    for filerep in filerep_lst:
    #        for tid in include_set:
    #            abun_cnt[tid] += weight*abund_dic[filerep][tid]
    output.write('\033[92m OK! \033[0m\n')
    if abundances:
        output.write(f'  \033[90mIncluding {len(abundances)} taxa at '
                     f'taxlevel "{level.value}"\033[0m\n'
                     '  \033[90mBuilding taxonomy tree...\033[0m')
    tree = TaxTree()
    tree.grow(children, _ROOT, [], abundances, {})
    output.write('\033[92m OK! \033[0m\n')
    # Prune the tree
    output.write('  \033[90mPruning taxonomy tree...\033[0m')
    prune_tree(tree, nmin=mintaxa, collapse=kollapse)
    output.write('\033[92m OK! \033[0m\n')
    # Get the taxa with their abundances and taxonomical levels
    output.write('  \033[90mFiltering taxa...\033[0m')
    abundances = col.Counter()
    taxlevels = {}
    get_taxa(tree, abundances, taxlevels,
             mindepth=0, maxdepth=0,
             include=including, exclude=excluding)
    abundances = +abundances  # remove zero and negative counts
    output.write('\033[92m OK! \033[0m\n')
    # Save the taxonomy in tsv file
    filen_shared = 'SHARED_' + level.value + _LINEAGEXT
    log = write_lineage(parents, names, tree,
                        filen_shared, abundances, kollapse)
    output.write(log)
    filen_lst.append(filen_shared)
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
        help='increase output verbosity (from -v to -vvv)'
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
            pool.close()
            map(mp.pool.ApplyResult.wait, async_results)
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
                    TaxLevel.selected_taxlevels))) as pool:
                async_results = [pool.apply_async(
                    process_taxlevel,
                    args=[level],
                    kwds=kwargs
                ) for level in TaxLevel.selected_taxlevels]
                pool.close()
                map(mp.pool.ApplyResult.wait, async_results)
                for level, (filen_diff_dic[level]) in zip(
                        TaxLevel.selected_taxlevels,
                        [r.get() for r in async_results]):
                    pass
        else:
            for level in TaxLevel.selected_taxlevels:
                filen_diff_dic[level] = process_taxlevel(level, **kwargs)

    # Generate the Krona html file calling ktImportText
    subprc = ["ktImportText"]
    subprc.extend(filen_lst)
    try:
        subprc.extend([filen_diff_dic[level][i]
                       for level in TaxLevel.selected_taxlevels
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
    main()

# Other useful functions or examples
#    print_tree(tree)
#    print_nodes(tree, nmin=1)
#    print(nodecnt.most_common(10))
#    nodeset = set(root)
#    get_taxa(tree, nodeset, maxdepth=3)
#    write_tsv(filerep + '.tsv', nodecnt)
#
# def get_taxa_bkp(subtree, node_container, mindepth=0, maxdepth=0, branch=False):
#     """
#     Recursive function to get the taxa between min and max depth levels
#
#     mindepth 0 (default) gets the taxa from the very beginning depth level.
#     maxdepth 0 (default) does not stop the search at any depth level.
#     node_container is a input/output set or col.Counter: at entry it should contain
#         the root of the subtree, while at output it will contain all the nodes
#         from the root to the maxdepth level (if any).
#     branch controls whether the root of the subtree in node_container at entry
#         should be considered the root of the branch whose taxa is to be got or
#         not (default).
#     """
#     _isset = False
#     if isinstance(node_container, set):
#         _isset = True
#     elif isinstance(node_container, col.Counter):
#         pass
#     else:
#         raise Exception('\n\033[91mERROR!\033[0m get_taxa: ' +
#                         'container for taxa not supported (use set or col.Counter)')
#     mindepth -= 1
#     maxdepth -= 1
#     for tid in subtree:
#         if subtree[tid]:
#             if (not branch) or (tid in node_container):
#                 if not maxdepth:
#                     break
#                 if mindepth <= 0:
#                     for child in subtree[tid]:
#                         if _isset:
#                             node_container.add(child)
#                         else:
#                             node_container[child] = subtree[tid][child].counts
#             get_taxa(subtree[tid], node_container, mindepth, maxdepth, branch)
#
# def write_tsv(tsvfile, node_cnt):
#     """Write tab separated values file with taxId and abundances."""
#     print('\033[90mSaving tsv file', tsvfile, '...\033[0m', end='')
#     sys.stdout.flush()
#     with open(tsvfile, 'w', newline='') as tsv_handle:
#         tsvwriter = csv.writer(tsv_handle, dialect='krona')
#         tsvwriter.writerow(["#taxID", "Abundance"])
#         for tax_id in node_cnt:
#             tsvwriter.writerow([tax_id, node_cnt[tax_id]])
#     print('\033[92m OK! \033[0m')
