"""
Taxonomy class, currently representing the Gene Ontology.

"""
import sys
from typing import Set, Counter, Iterable, Tuple

from recentrifuge.config import Filename, Id, Parents, Names, Children
from recentrifuge.config import GO_ROOT, CELLULAR_ORGANISMS
from recentrifuge.config import gray, green
from recentrifuge.ontology import Ontology
from recentrifuge.rank import Ranks, Rank, UnsupportedTaxLevelError
from recentrifuge.trees import TaxTree

GOid = Id


class GeneOntology(Ontology):
    """Taxonomy related data and methods."""

    def __init__(self,
                 nodes_file: Filename,
                 names_file: Filename,
                 excluding: Set[Id] = None,
                 including: Set[Id] = None,
                 debug: bool = False,
                 ) -> None:

        # Type data declaration and initialization
        self.ROOT = GO_ROOT
        self.parents = Parents({})
        self.children = Children({})
        self.ranks: Ranks = Ranks({})
        self.names: Names = Names({})
        self.debug: bool = debug

        # Initialization methods
        self.read_nodes(nodes_file)
        self.read_names(names_file)
        self.build_children()
        self.build_ranks()

        # Show explicitly included and excluded taxa
        if including:
            print('List of GOs (and below) to be explicitly included:')
            print('\t\tGO id\tName')
            for goid in including:
                print(f'\t\t{goid}\t{self.names[goid]}')
        else:
            # To excluding to operate not on single taxa but on subtrees
            including = {GO_ROOT}
        self.including: Set[Id] = including
        if excluding:
            print('List of GOs (and below) to be excluded:')
            print('\t\tGO id\tName')
            for goid in excluding:
                print(f'\t\t{goid}\t{self.names[goid]}')
        self.excluding: Set[Id] = excluding

    def read_nodes(self, nodes_file: Filename) -> None:
        """Build dicts of parent and rank for a given GOid (key)"""
        print('\033[90mLoading GO nodes...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(nodes_file, 'r') as file:
                file.readline()  # Discard header
                for line in file:
                    _goid, _parent = line.rstrip().split('\t')
                    goid = GOid(_goid)
                    parent = GOid(_parent)
                    self.parents[goid] = parent
        except OSError:
            print(f'\n\033[91mERROR!\033[0m Cannot read {nodes_file}')
            raise
        except ValueError:
            print(f'\n\033[91mERROR!\033[0m Problem with {nodes_file}, line:')
            print(line)
            raise
        else:
            # Create relations to GO_ROOT
            self.parents[GOid('GO:0003674')] = GO_ROOT  # Molecular function
            self.parents[GOid('GO:0005554')] = GO_ROOT  # ... alternate ID
            self.parents[GOid('GO:0005575')] = GO_ROOT  # Cellular component
            self.parents[GOid('GO:0008372')] = GO_ROOT  # ... alternate ID
            self.parents[GOid('GO:0008150')] = GO_ROOT  # Biological process
            self.parents[GOid('GO:0000004')] = GO_ROOT  # ... alternate ID
            self.parents[GOid('GO:0007582')] = GO_ROOT  # ... alternate ID
            self.parents[GOid('GO:0044699')] = GO_ROOT  # ... alternate ID
            print('\033[92m OK! \033[0m')

    def read_names(self, names_file: Filename) -> None:
        """Build dict with name for a given GOid (key)."""
        print('\033[90mLoading GO names...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(names_file, 'r') as file:
                file.readline()  # Discard header
                for line in file:
                    _goid, _activity, _name = line.rstrip().split('\t')
                    goid = Id(_goid)
                    self.names[goid] = _name
                    self.ranks[goid] = Rank.NO_RANK
        except OSError:
            raise Exception('\n\033[91mERROR!\033[0m Cannot read "' +
                            names_file + '"')
        else:
            self.names[GO_ROOT] = str(GO_ROOT)
            print('\033[92m OK! \033[0m')

    def build_children(self) -> None:
        """Build dict of children for a given parent GOid (key)."""
        print('\033[90mBuilding dict of parent to children nodes...\033[0m',
              end='')
        sys.stdout.flush()
        self.children: Children = Children({})
        for goid in self.parents:
            if self.parents[goid] not in self.children:
                self.children[self.parents[goid]] = {}
            self.children[self.parents[goid]][goid] = 0
        print('\033[92m OK! \033[0m')

    def build_ranks(self) -> None:
        """Build virtual rank (vrank) hierarchy for gene ontology. """
        print(gray('Building virtual rank hierarchy for gene ontology...'),
              end='')
        tree = TaxTree()
        # Grow tree and shape it
        tree.grow(ontology=self, look_ancestors=False)
        tree.vrank()
        # Get rank data from tree
        tree.get_taxa(ranks=self.ranks)
        print(green(' OK!'))

    def get_rank(self, goid: Id) -> Rank:
        """Retrieve the rank for a GOid."""
        return self.ranks.get(goid, Rank.NO_RANK)

    def get_name(self, goid: Id) -> str:
        """Retrieve the name for a GOid."""
        return self.names.get(goid, 'Unnamed')

    def get_ancestors(self, leaves: Iterable[GOid]
                      ) -> Tuple[Set[GOid], Set[GOid]]:
        """Return the GOids entered with all their ancestors"""
        ancestors: Set[GOid] = set(leaves)
        orphans: Set[GOid] = set()
        for leaf in leaves:
            goid: GOid = leaf
            while goid != GO_ROOT:
                try:
                    goid = self.parents[goid]
                except KeyError:
                    orphans.add(goid)
                    break
                else:
                    ancestors.add(goid)
        return ancestors, orphans
