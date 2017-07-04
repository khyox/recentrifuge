"""
Recentrifuge core classes and functions.

"""
# pylint: disable=not-an-iterable, bad-reversed-sequence, E0203
import collections as col
import copy
import csv
import io
import subprocess
import sys
from enum import Enum
from typing import List, Iterator, Set, Counter, Iterable, Tuple, Union, Dict
from typing import NewType

from recentrifuge.config import Filename, Sample, TaxId, Children, Names
from recentrifuge.config import Parents, ROOT, HTML_SUFFIX, CELLULAR_ORGANISMS
from recentrifuge.krona import COUNT, UNASSIGNED, TID, RANK, SCORE
from recentrifuge.krona import KronaTree, Elm

# Type annotations
# pylint: disable=invalid-name
# Ranks and Levels are devised to be one the inverse of the other
Ranks = NewType('Ranks', Dict[TaxId, 'Rank'])  # Rank of each TaxId
TaxLevels = NewType('TaxLevels', Dict['Rank', Set[TaxId]])  # TaxIds for rank


# pylint: enable=invalid-name

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
        obj._value_ = _value  # pylint: disable=protected-access
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
                 collapse: bool = True,
                 excluding: Set[TaxId] = None,
                 including: Set[TaxId] = None,
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
            including = {ROOT}
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
                    if self.collapse and parent == CELLULAR_ORGANISMS:
                        self.parents[tid] = ROOT
                    else:
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

    def get_name(self, taxid: TaxId) -> str:
        """Retrieve the name for a TaxId."""
        return self.names.get(taxid, 'Unnamed')


class TaxTree(dict):
    """Nodes of a taxonomical tree"""

    def __init__(self, *args,
                 counts: int = 0,
                 rank: Rank = Rank.UNCLASSIFIED,
                 score: float = 0.5
                 ) -> None:
        super().__init__(args)
        self.counts: int = counts
        self.taxlevel: Rank = rank
        self.score: float = score
        # Accumulated counts are optionally populated with self.accumulate()
        self.acc: int = 0

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
             taxid: TaxId = ROOT,
             _path: List[TaxId] = None,
             ) -> None:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            abundance: counter for taxids with their abundances.
            taxid: It's ROOT by default for the first/base method call
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
                    self[taxid].grow(taxonomy=taxonomy,
                                     abundance=abundance,
                                     taxid=child,
                                     _path=_path + [taxid])

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
                if target in self[tid] and target != ROOT:
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
            if tid == ROOT:
                nodes_traced[ROOT] = [ROOT, ]  # Root node special case
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
                 abundance: Counter[TaxId] = None,
                 accs: Counter[TaxId] = None,
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

        The three outputs are optional, abundance, accumulated
        abundance and rank: to enable them, an empty dictionary should
        be attached to them. Makes no sense to call the method without,
        at least, one of them.

        Args:
            abundance: Optional I/O dict; at 1st entry should be empty.
            accs: Optional I/O dict; at 1st entry should be empty.
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
                    if abundance is not None:
                        abundance[tid] = self[tid].counts
                    if accs is not None:
                        accs[tid] = self[tid].acc
                    if ranks is not None:
                        ranks[tid] = self[tid].taxlevel
                if self[tid]:
                    self[tid].get_taxa(abundance, accs, ranks,
                                       mindepth, maxdepth,
                                       include, exclude,
                                       just_level, in_branch)

    def toxml(self,
              taxonomy: Taxonomy,
              krona: KronaTree,
              node: Elm = None,
              mindepth: int = 0,
              maxdepth: int = 0,
              include: Union[Tuple, Set[TaxId]] = (),
              exclude: Union[Tuple, Set[TaxId]] = (),
              _in_branch: bool = False
              ) -> None:
        """
        Recursively convert to XML between min and max depth levels.

        Args:
            taxonomy: Taxonomy object.
            krona: Input/Output KronaTree object to be populated.
            node: Base node (None to use the root of krona argument).
            mindepth: 0 gets the taxa from the very beginning depth.
            maxdepth: 0 does not stop the search at any depth level.
            include: contains the root taxid of the subtrees to be
                included. If it is empty (default) all the taxa is
                included (except explicitly excluded).
            exclude: contains the root taxid of the subtrees to be
                excluded
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
                new_node: Elm
                if mindepth <= 0 and in_branch:
                    if node is None:
                        node = krona.getroot()
                    new_node = krona.node(
                        node, taxonomy.get_name(tid),
                        {COUNT: {krona.samples[0]: str(self[tid].acc)},
                         UNASSIGNED: {krona.samples[0]: str(self[tid].counts)},
                         TID: str(tid),
                         RANK: taxonomy.get_rank(tid).name.lower(),
                         SCORE: {krona.samples[0]: str(self[tid].score)}}
                    )
                if self[tid]:
                    self[tid].toxml(taxonomy,
                                    krona, new_node,
                                    mindepth, maxdepth,
                                    include, exclude,
                                    in_branch)

    def accumulate(self) -> None:
        """
        Recursively populate accumulated counts (self.acc).

        Accumulate counts in higher taxonomical levels.
        """
        self.acc = self.counts  # Initialize accumulated with unassigned counts
        if self:  # If this node has subtrees
            for tid in list(self):
                self[tid].accumulate()  # Accumulate for each branch/leaf
                self.acc += self[tid].acc  # Acc lower tax acc in this node

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


class MultiTree(dict):
    """Nodes of a multiple taxonomical tree"""

    def __init__(self, *args,
                 samples: List[Sample],
                 counts: Dict[Sample, int] = None,
                 accs: Dict[Sample, int] = None,
                 rank: Rank = Rank.UNCLASSIFIED,
                 scores: Dict[Sample, float] = None
                 ) -> None:
        """

        Args:
            *args: Arguments to parent class (dict)
            samples: List of samples to coexist in the (XML) tree
            counts: Dict of abundance for each sample in this tax node
            accs: Dict of accumulated abundance for each sample
            rank: Rank of this tax node
            scores: Dict with score for each sample
        """
        super().__init__(args)
        self.samples: List[Sample] = samples
        self.taxlevel: Rank = rank
        # Dict(s) to List(s) following samples order to save space in each node
        if counts is None:
            counts = {sample: 0 for sample in samples}
        self.counts: List[int] = [counts[sample] for sample in samples]
        if accs is None:
            accs = {sample: 0 for sample in samples}
        self.accs: List[int] = [accs[sample] for sample in samples]
        if scores is None:
            scores = {sample: 1.0 for sample in samples}
        self.score: List[float] = [scores[sample] for sample in samples]

    def __str__(self, num_min: int = 1) -> None:
        """
        Recursively print populated nodes of the taxonomy tree

        Args:
            num_min: minimum abundance of a node to be printed

        Returns: None

        """
        for tid in self:
            if max(self[tid].counts) >= num_min:
                print(f'{tid}[{self[tid].counts}]', end='')
            if self[tid]:  # Has descendants
                print('', end='->(')
                self[tid].__str__(num_min=num_min)
            else:
                print('', end=',')
        print(')', end='')

    def grow(self,
             taxonomy: Taxonomy,
             abundances: Dict[Sample, Counter[TaxId]],
             accs: Dict[Sample, Counter[TaxId]],
             taxid: TaxId = ROOT,
             _path: List[TaxId] = None) -> None:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            abundances: Dict of counters with taxids' abundances.
            accs: Dict of counters with taxids' accumulated abundances.
            taxid: It's ROOT by default for the first/base method call
            _path: list used by the recursive algorithm to avoid loops

        Returns: None

        """
        if not _path:
            _path = []
        if taxid not in _path:  # Avoid loops for repeated taxid (like root)
            multicount: Dict[Sample, int] = {
                sample: abundances[sample].get(taxid, 0)
                for sample in self.samples
            }
            multiacc: Dict[Sample, int] = {
                sample: accs[sample].get(taxid, 0)
                for sample in self.samples
            }
            if any(multiacc.values()):  # Check for any populated branch
                self[taxid] = MultiTree(samples=self.samples,
                                        counts=multicount,
                                        accs=multiacc,
                                        rank=taxonomy.get_rank(taxid))
                if taxid in taxonomy.children:  # taxid has children
                    for child in taxonomy.children[taxid]:
                        self[taxid].grow(taxonomy=taxonomy,
                                         abundances=abundances,
                                         accs=accs,
                                         taxid=child,
                                         _path=_path + [taxid])

    def toxml(self,
              taxonomy: Taxonomy,
              krona: KronaTree,
              node: Elm = None,
              ) -> None:
        """
        Recursively get the taxa between min and max depth levels.

        Args:
            taxonomy: Taxonomy object.
            krona: Input/Output KronaTree object to be populated.
            node: Base node (None to use the root of krona argument).

        Returns: None

        """
        for tid in self:
            if node is None:
                node = krona.getroot()
            num_samples = len(self.samples)
            new_node: Elm = krona.node(
                parent=node,
                name=taxonomy.get_name(tid),
                values={COUNT: {self.samples[i]: str(self[tid].accs[i])
                                for i in range(num_samples)},
                        UNASSIGNED: {self.samples[i]: str(self[tid].counts[i])
                                     for i in range(num_samples)},
                        TID: str(tid),
                        RANK: taxonomy.get_rank(tid).name.lower(),
                        SCORE: {self.samples[i]: str(self[tid].score[i])
                                for i in range(num_samples)},
                        }
            )
            if self[tid]:
                self[tid].toxml(taxonomy=taxonomy,
                                krona=krona,
                                node=new_node)


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

def process_rank(*args, **kwargs):
    """
    Process results for a taxlevel (to be usually called in parallel!).
    """

    # Recover input and parameters
    rank: Rank = args[0]
    taxonomy: Taxonomy = kwargs['taxonomy']
    including = taxonomy.including
    excluding = taxonomy.excluding
    mintaxa = kwargs['mintaxa']
    trees: Dict[Sample, TaxTree] = kwargs['trees']
    taxids: Dict[Sample, TaxLevels] = kwargs['taxids']
    reports: List[Sample] = kwargs['reports']
    control: bool = kwargs['control']

    # Declare/define variables
    samples: List[Sample] = []
    sample: Sample
    abundances: Dict[Sample, Counter[TaxId]] = {}
    accumulators: Dict[Sample, Counter[TaxId]] = {}
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
        leveled_tree.accumulate()
        # Get abundance for the exclusive analysis
        exclusive_abundance: Counter[TaxId] = col.Counter()
        exclusive_acc: Counter[TaxId] = col.Counter()
        leveled_tree.get_taxa(exclusive_abundance, exclusive_acc,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=exclude,  # Extended exclusion set
                              just_level=rank,
                              )
        # Generate a new tree to be able to calculate accs
        tree = TaxTree()
        tree.grow(taxonomy, exclusive_abundance)
        tree.accumulate()  # Calculate the accumulated values
        exclusive_acc = col.Counter()
        tree.get_taxa(None, exclusive_acc,
                      include=including,
                      exclude=excluding)
        output.write('\033[92m OK! \033[0m\n')
        sample = Sample(f'EXCLUSIVE_{rank.name.lower()}_{report}')
        samples.append(sample)
        abundances[sample] = exclusive_abundance
        accumulators[sample] = exclusive_acc
        # Get partial abundance for the shared analysis
        sub_shared_abundance: SharedCounter = SharedCounter()
        leveled_tree.get_taxa(sub_shared_abundance,
                              mindepth=0, maxdepth=0,
                              include=including,
                              exclude=excluding,
                              just_level=rank,
                              )
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
        tree.accumulate()  # Calculate the accumulated values
        new_abund: Counter[TaxId] = col.Counter()
        new_accs: Counter[TaxId] = col.Counter()
        tree.get_taxa(new_abund, new_accs,
                      mindepth=0, maxdepth=0,
                      include=including,
                      exclude=excluding)
        output.write('\033[92m OK! \033[0m\n')
        sample = Sample(f'SHARED_{rank.name.lower()}')
        samples.append(sample)
        abundances[Sample(sample)] = new_abund
        accumulators[Sample(sample)] = new_accs

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
            leveled_tree.accumulate()  # Calculate the accumulated values
            control_abundance: Counter[TaxId] = col.Counter()
            leveled_tree.get_taxa(control_abundance,
                                  include=including,
                                  exclude=exclude,  # Extended exclusion set
                                  just_level=rank,
                                  )
            # Generate a new tree to be able to calculate accs
            tree = TaxTree()
            tree.grow(taxonomy, control_abundance)
            tree.accumulate()  # Calculate the accumulated values
            control_acc: Counter[TaxId] = col.Counter()
            tree.get_taxa(None, control_acc,
                          include=including,
                          exclude=excluding)
            output.write('\033[92m OK! \033[0m\n')
            control_abundance = +control_abundance  # remove counts <= 0
            sample = Sample(f'CONTROL_{rank.name.lower()}_{report}')
            samples.append(sample)
            abundances[sample] = control_abundance
            accumulators[sample] = control_acc

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
            tree.accumulate()
            new_shared_ctrl_abundance: Counter[TaxId] = col.Counter()
            new_shared_ctrl_acc: Counter[TaxId] = col.Counter()
            tree.get_taxa(new_shared_ctrl_abundance,
                          include=including,
                          exclude=exclude,  # Extended exclusion set
                          just_level=rank,
                          )
            # Generate a new tree to be able to calculate accs
            tree = TaxTree()
            tree.grow(taxonomy, new_shared_ctrl_abundance)
            tree.accumulate()  # Calculate the accumulated values
            tree.get_taxa(None, new_shared_ctrl_acc,
                          include=including,
                          exclude=excluding)

            output.write('\033[92m OK! \033[0m\n')
            sample = Sample(f'SHARED_CONTROL_{rank.name.lower()}')
            samples.append(sample)
            abundances[sample] = new_shared_ctrl_abundance
            accumulators[sample] = new_shared_ctrl_acc

    # Print output and return
    print(output.getvalue())
    sys.stdout.flush()
    return samples, abundances, accumulators


def write_lineage(parents: Parents,
                  names: Dict[TaxId, str],
                  tree: TaxTree,
                  lineage_file: str,
                  nodes: Counter[TaxId],
                  collapse: bool = True,
                  ) -> str:
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
                    taxids_dic[tid].remove(CELLULAR_ORGANISMS)
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


def krona_from_text(samples: List[Sample],
                    outputs: Dict[Rank, Filename],
                    htmlfile: Filename = Filename('Output' + HTML_SUFFIX),
                    ):
    """Generate the Krona html file calling ktImportText.

    Superseded by krona.krona_from_xml().

    """
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


