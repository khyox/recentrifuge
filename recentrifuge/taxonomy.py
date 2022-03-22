"""
Taxonomy class, currently representing the NCBI taxonomy.

"""
import collections as col
import re
import sys
from typing import Set, Counter, Iterable, Tuple, Union, List, Dict

from recentrifuge.config import Filename, Id, Parents, Names, Children
from recentrifuge.config import ROOT, CELLULAR_ORGANISMS, UNNAMED
from recentrifuge.config import red, magenta
from recentrifuge.rank import Ranks, Rank, UnsupportedTaxLevelError
from recentrifuge.ontology import Ontology


class Taxonomy(Ontology):
    """Taxonomy related data and methods."""

    # pylint: disable=unsubscriptable-object,unsupported-assignment-operation
    # pylint: disable=not-an-iterable,unsupported-membership-test,no-member

    def __init__(self,
                 nodes_file: Filename,
                 names_file: Filename,
                 plasmid_file: Filename = None,
                 collapse: bool = True,
                 excluding: Union[Tuple, Set[Id]] = (),
                 including: Union[Tuple, Set[Id]] = (),
                 debug: bool = False,
                 ) -> None:

        # Type data declaration and initialization
        self.ROOT = ROOT
        self.parents: Parents = Parents({})
        self.ranks: Ranks = Ranks({})
        self.names: Names = Names({})
        self.children: Children = Children({})
        self.collapse: bool = collapse
        self.debug: bool = debug

        # Initialization methods
        self.read_nodes(nodes_file)
        self.read_names(names_file)
        if plasmid_file:
            self.read_plasmids(plasmid_file)
        self.build_children()

        # Show explicitly included and excluded taxa
        if including:
            print('List of taxa (and below) to be explicitly included:')
            print('\t\tId\tScientific Name')
            for taxid in including:
                try:
                    print(f'\t\t{taxid}\t{self.names[taxid]}')
                except KeyError:
                    print(red('ERROR!'), 'The taxid', taxid,
                          'to include seems invalid.', red('Aborting!'))
                    raise
        else:
            # For excluding to operate not on single taxa but on subtrees
            including = {ROOT}
        self.including: Union[Tuple, Set[Id]] = including
        if excluding:
            print('List of taxa (and below) to be excluded:')
            print('\t\tId\tScientific Name')
            for taxid in excluding:
                try:
                    print(f'\t\t{taxid}\t{self.names[taxid]}')
                except KeyError:
                    print(red('ERROR!'), 'The taxid', taxid,
                          'to exclude seems invalid.', red('Aborting!'))
                    raise
        self.excluding: Union[Tuple, Set[Id]] = excluding

    def read_nodes(self, nodes_file: Filename) -> None:
        """Build dicts of parent and rank for a given taxid (key)"""
        print('\033[90mLoading NCBI nodes...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(nodes_file, 'r') as file:
                for line in file:
                    _tid, _parent, _rank, *_ = line.split('\t|\t')
                    tid = Id(_tid)
                    parent = Id(_parent)
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
            print(red('ERROR!'), f'Cannot read {nodes_file}.')
            print(magenta('TIP:'),
                  'Did you select the right path with the "-n" option?')
            print(magenta('TIP:'),
                  'Did you use "Retaxdump" to install the dump files?')
            raise
        else:
            print('\033[92m OK! \033[0m')

    def read_names(self, names_file: Filename) -> None:
        """Build dict with name for a given taxid (key)."""
        print('\033[90mLoading NCBI names...\033[0m', end='')
        sys.stdout.flush()
        try:
            with open(names_file, 'r') as file:
                for line in file:
                    if 'scientific name' in line:
                        tid, scientific_name, *_ = line.split('\t|\t')
                        self.names[Id(tid)] = scientific_name
        except OSError:
            print(red('ERROR!'), f'Cannot read {names_file}.')
            print(magenta('TIP:'),
                  'Did you use "Retaxdump" to install the dump files?')
            raise
        else:
            print('\033[92m OK! \033[0m')

    def read_plasmids(self, plasmid_file: Filename) -> None:
        """Read, check and include plasmid data"""
        print('\033[90mLoading LMAT plasmids...\033[0m', end='')
        sys.stdout.flush()
        pattern1 = re.compile(
            r"""((?:"([\w\-\.\(\)/+=':,%\*\s]*)"$)|(?:^([\w\-\.\(\)/+=':\*\s]*(?:, (?:strain|isolate|plasmid) [\w\-/\.]*)*(?:, fragment \w*)?(?:, contig \w)?)(?:, a cloning vector)?(?=(?=(?:, (?:complete|partial) (?:plasmid |genomic )*(?:sequence|genome|cds|replicon))(?:\[sequence_id)*)|(?:, complete sequence)*, whole genome shotgun sequence|\[sequence_id)))"""  # pylint: disable=line-too-long
        )
        pattern2 = re.compile(
            r"""(^(?:[A-Za-z0-9/=\-\.{},]*(?: |.)){1,8})"""
        )
        match: Counter = col.Counter()
        try:
            with open(plasmid_file, 'r') as file:
                for line in file:
                    _tid, _parent, *_, last = line.rstrip('\n').split('\t')
                    last = last.split(r'|')[-1]
                    tid = Id(_tid)
                    parent = Id(_parent)
                    # Plasmids sanity checks
                    if tid in self.parents:  # if plasmid tid already in NCBI
                        match['ERR1'] += 1
                        if self.debug:
                            print(f'\033[93mPlasmid taxid ERROR!\033[0m'
                                  f' Taxid={tid} already a NCBI taxid. '
                                  f'Declared parent is {parent} but '
                                  f'NCBI parent is {self.parents[tid]}.')
                            print('\tPlasmid details: ', last)
                        continue
                    elif tid == parent:  # if plasmid and parent tids are equal
                        match['ERR2'] += 1
                        if self.debug:
                            print(f'\033[93mPlasmid parent taxid ERROR!\033[0m'
                                  f' Taxid={tid} and parent={parent}.')
                            print('\t\t   Plasmid details: ', last)
                        continue
                    else:  # No problem, go ahead and add the plasmid!
                        self.parents[tid] = parent
                    # Plasmid name extraction by regular expressions
                    name: str
                    try:
                        name = pattern1.search(last).group(1)  # type: ignore
                        name = 'Plasmid ' + name.strip(r'"').strip(',')
                    except AttributeError:
                        try:
                            name = pattern2.search(  # type: ignore
                                last).group(1).strip()
                            name = 'Plasmid ' + name
                        except AttributeError:
                            name = 'Plasmid ' + tid
                            match['FAIL'] += 1
                        else:
                            match['PAT2'] += 1
                    else:
                        match['PAT1'] += 1
                    self.names[tid] = name
        except OSError:
            print('\033[93mWARNING\033[0m: Cannot read "' +
                        plasmid_file + '". Plasmid taxids not loaded!')
            print(magenta('TIP:'),
                  'Manual installation of the plasmids file required.')
            raise
        else:  # Statistics about plasmids
            print('\033[92m OK! \033[0m\n',
                  '\033[90mPlasmid sanity check:\033[0m',
                  f'\033[93m rejected\033[0m (taxid error) = {match["ERR1"]}',
                  f'\033[93m rejected\033[0m (parent error) = {match["ERR2"]}')
            print('\033[90m Plasmid pattern matching:\033[0m',
                  f'\033[90m 1st type =\033[0m {match["PAT1"]} ',
                  f'\033[90m 2nd type =\033[0m {match["PAT2"]} ',
                  f'\033[90m other =\033[0m {match["FAIL"]}')

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

    def get_rank(self, taxid: Id) -> Rank:
        """Retrieve the rank for a Id."""
        return self.ranks.get(taxid, Rank.UNCLASSIFIED)

    def get_name(self, taxid: Id) -> str:
        """Retrieve the name for a Id."""
        return self.names.get(taxid, UNNAMED)

    def get_ancestors(self, leaves: Iterable[Id]
                      ) -> Tuple[Set[Id], Set[Id]]:
        """Return the taxids entered with all their ancestors"""
        ancestors: Set[Id] = set(leaves)
        orphans: Set[Id] = set()
        for leaf in leaves:
            tid: Id = leaf
            while tid != ROOT:
                try:
                    tid = self.parents[tid]
                except KeyError:
                    orphans.add(tid)
                    break
                else:
                    ancestors.add(tid)
        return ancestors, orphans

    def get_lineage(self, taxid: Id, as_dict: bool = False,
                    add_root: bool = False) -> List[Tuple[Rank, Id, str]]:
        """Retrieve the lineage of a taxon up to the root of the taxonomy"""
        lineage: List[Tuple[Rank, Id, str]] = []
        lineage_as_dict: Dict[Rank, Dict[str, str]] = {}
        if (self.get_rank(taxid) is Rank.UNCLASSIFIED
                or self.get_name(taxid) is UNNAMED):
            print(red('ERROR!'), f' Taxid {taxid} seems not valid.')
            raise KeyError('Invalid taxonomic id')
        no_rank_cnt: int = 1
        while taxid != ROOT:
            taxid_label: str = 'taxid'
            name_label: str = 'name'
            this_rank: Rank = self.get_rank(taxid)
            if this_rank is Rank.NO_RANK:
                taxid_label = 'taxid_' + str(no_rank_cnt)
                name_label = 'name_' + str(no_rank_cnt)
                no_rank_cnt += 1
            lineage.append((this_rank, taxid, self.get_name(taxid)))
            if lineage_as_dict.get(this_rank) is None:
                lineage_as_dict[this_rank] = {}
            lineage_as_dict[this_rank].update({
                taxid_label: taxid, name_label: self.get_name(taxid)})
            try:
                taxid = self.parents[taxid]
            except KeyError:
                print(red('ERROR!'), f'Orphan taxid {taxid} found.')
                raise
        if add_root:
            root_rank: Rank = self.get_rank(ROOT)
            lineage.append((root_rank, ROOT, self.get_name(ROOT)))
            if lineage_as_dict.get(root_rank) is None:
                lineage_as_dict[root_rank] = {}
            lineage_as_dict[root_rank].update({
                'taxid': ROOT, 'name': self.get_name(ROOT)})
        if as_dict:
            return lineage_as_dict
        else:
            return lineage
