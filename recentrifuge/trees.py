"""
TaxTree and MultiTree classes.

"""

import collections as col
import io
from typing import Counter, Union, Dict, List, Iterable, Tuple, Set
from recentrifuge.config import TaxId, Parents, Sample, Score, Scores
from recentrifuge.config import ROOT, NO_SCORE, UnionCounter, UnionScores
from recentrifuge.krona import COUNT, UNASSIGNED, TID, RANK, SCORE
from recentrifuge.krona import KronaTree, Elm
from recentrifuge.rank import Rank, Ranks, TaxLevels
from recentrifuge.shared_counter import SharedCounter
from recentrifuge.taxonomy import Taxonomy


class SampleDataByTaxId(object):
    """Typical data in a sample ordered by taxonomical id"""

    def __init__(self, init: List[str] = None) -> None:
        """Initialize data structures

        Individual options: 'counts', 'ranks', 'scores', 'accs',
          'shared_counts', 'shared_scores'
        Group options: 'all' initialize all non-shared, and
           'shared' initializes both counts and scores with
            SharedCounter."""
        self.counts: UnionCounter = None
        self.ranks: Ranks = None
        self.scores: UnionScores = None
        self.accs: Counter[TaxId] = None
        if 'counts' in init or 'all' in init:
            self.counts = Counter()
        if 'ranks' in init or 'all' in init:
            self.ranks = Ranks({})
        if 'scores' in init or 'all' in init:
            self.scores = Scores({})
        if 'accs' in init or 'all' in init:
            self.accs = Counter()
        if 'shared_counts' in init or 'shared' in init:
            self.counts = SharedCounter()
        if 'shared_scores' in init or 'shared' in init:
            self.scores = SharedCounter()

    def set(self,
            counts: UnionCounter = None,
            ranks: Ranks = None,
            scores: UnionScores = None,
            accs: Counter[TaxId] = None) -> None:
        """Set the data fields"""
        if counts is not None:
            self.counts = counts
        if ranks is not None:
            self.ranks = ranks
        if scores is not None:
            self.scores = scores
        if accs is not None:
            self.accs = accs

    def get_counts(self) -> Counter[TaxId]:
        """Get (non shared) counts"""
        if isinstance(self.counts, Counter):
            return self.counts
        raise TypeError

    def get_shared_counts(self) -> SharedCounter:
        """Get shared counts"""
        if isinstance(self.counts, SharedCounter):
            return self.counts
        raise TypeError

    def get_scores(self) -> Scores:
        """Get (non shared) scores"""
        if isinstance(self.scores, dict):
            return self.scores  # type: ignore
        raise TypeError

    def get_shared_scores(self) -> SharedCounter:
        """Get shared scores"""
        if isinstance(self.scores, SharedCounter):
            return self.scores
        raise TypeError

    def get_accs(self) -> Counter[TaxId]:
        """Get accumulated counter"""
        if isinstance(self.accs, Counter):
            return self.accs
        raise TypeError

    def clear(self, fields: List[str] = None) -> None:
        """Clear the data field"""
        if 'counts' in fields or 'all' in fields and self.counts is not None:
            self.counts.clear()
        if 'ranks' in fields or 'all' in fields and self.ranks is not None:
            self.ranks.clear()
        if 'scores' in fields or 'all' in fields and self.scores is not None:
            self.scores.clear()
        if 'accs' in fields or 'all' in fields and self.accs is not None:
            self.accs.clear()

    def purge_counters(self) -> None:
        """Purge elements with zero counts in counters"""
        if isinstance(self.counts, Counter):
            self.counts = +self.counts  # pylint: disable=E1130
        if isinstance(self.accs, Counter):
            self.accs = +self.accs  # pylint: disable=E1130

    def get_taxlevels(self) -> TaxLevels:
        """Get TaxLevels (taxids of ranks) from Ranks (rank of taxids)"""
        if self.ranks:
            return Rank.ranks_to_taxlevels(self.ranks)
        return NotImplemented


class TaxTree(dict):
    """Nodes of a taxonomical tree"""

    def __init__(self, *args,
                 counts: int = 0,
                 rank: Rank = Rank.UNCLASSIFIED,
                 score: float = 0,
                 acc: int = 0,
                 ) -> None:
        super().__init__(args)
        self.counts: int = counts
        self.rank: Rank = rank
        self.score: float = score
        self.acc: int = acc

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
             counts: Counter[TaxId] = None,
             scores: Union[Dict[TaxId, Score], SharedCounter] = None,
             ancestors: Set[TaxId] = None,
             taxid: TaxId = ROOT,
             _path: List[TaxId] = None,
             ) -> None:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            counts: counter for taxids with their abundances.
            scores: optional dict with the score for each taxid.
            ancestors: optional set of ancestors.
            taxid: It's ROOT by default for the first/base method call
            _path: list used by the recursive algorithm to avoid loops

        Returns: None

        """
        # Checks and initializations in the first call to the function
        if not _path:
            _path = []
            if not counts:
                abundances = col.Counter({ROOT: 1})
            if not scores:
                scores = {}
            if not ancestors:
                ancestors, _ = taxonomy.get_ancestors(counts.keys())

        # Go ahead if there is an ancestor or not repeated taxid (like root)
        if taxid in ancestors and taxid not in _path:
            self[taxid] = TaxTree(counts=counts.get(taxid, 0),
                                  score=scores.get(taxid, NO_SCORE),
                                  rank=taxonomy.get_rank(taxid))
            if taxid in taxonomy.children:  # taxid has children
                for child in taxonomy.children[taxid]:
                    self[taxid].grow(taxonomy=taxonomy,
                                     counts=counts,
                                     scores=scores,
                                     ancestors=ancestors,
                                     taxid=child,
                                     _path=_path + [taxid])

    def allin1(self,
               taxonomy: Taxonomy,
               counts: Counter[TaxId] = None,
               scores: Union[Dict[TaxId, Score], 'SharedCounter'] = None,
               ancestors: Set[TaxId] = None,
               tid: TaxId = ROOT,
               min_taxa: int = 1,
               min_rank: Rank = None,
               just_min_rank: bool = False,
               include: Union[Tuple, Set[TaxId]] = (),
               exclude: Union[Tuple, Set[TaxId]] = (),
               out: SampleDataByTaxId = None,
               _path: List[TaxId] = None) -> Union[int, None]:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            counts: counter for taxids with their abundances.
            scores: optional dict with the score for each taxid.
            ancestors: optional set of ancestors.
            tid: It's ROOT by default for the first/base method call
            min_taxa: minimum taxa to avoid pruning/collapsing
                one level to the parent one.
            min_rank: if any, minimum Rank allowed in the TaxTree.
            just_min_rank: If set, just min_rank taxa will be counted.
            include: contains the root taxid of the subtrees to be
                included. If it is empty (default) all the taxa is
                included (except explicitly excluded).
            exclude: root taxid of the subtrees to be excluded
            out: Optional I/O object, at 1st entry should be empty.
            _path: list used by the recursive algorithm to avoid loops

        Returns: Accumulated counts of new node (or None for no node)

        """

        def swmean(cnt1: int, sco1: Score, cnt2: int, sco2: Score) -> Score:
            """Weighted mean of scores by counts"""
            if sco1 == NO_SCORE:
                return sco2
            return Score((cnt1 * sco1 + cnt2 * sco2) / (cnt1 + cnt2))

        def populate_output(taxid: TaxId, source: TaxTree) -> None:
            """Populate the output structure"""
            nonlocal out
            if out.counts is not None:
                out.counts[taxid] = source[taxid].counts
            if out.ranks is not None:
                out.ranks[taxid] = source[taxid].rank
            if out.scores is not None and source[taxid].score != NO_SCORE:
                out.scores[taxid] = source[taxid].score
            if out.accs is not None:
                out.accs[taxid] = source[taxid].acc

        # Checks and initializations in the first call to the function
        if not _path:
            _path = []
            if not counts:
                counts = col.Counter({ROOT: 1})
            if not scores:
                scores = {}
            if min_rank is None and just_min_rank:
                raise RuntimeError('allin1: just_min_rank without min_rank')
            if not ancestors:
                ancestors, _ = taxonomy.get_ancestors(counts.keys())

        # Return if not an ancestor, or excluded taxa or loops for
        #  repeated taxid (like root); conditions ordered by prob
        if tid not in ancestors or tid in exclude or tid in _path:
            return None

        rank: Rank = taxonomy.get_rank(tid)
        parent_rank: Rank = None
        if min_rank:  # Get parent rank (NO_RANK is not an option)
            if tid == ROOT:
                parent_rank = Rank.ROOT
            else:
                for taxid in reversed(_path):
                    parent_rank = taxonomy.get_rank(taxid)
                    if parent_rank is not Rank.NO_RANK:
                        break
        abun: int = 0
        # For assigning counts, check just_min_rank conditions and,
        #  if include list, check if any include taxon is in _path
        if ((not just_min_rank or rank <= min_rank or parent_rank <= min_rank)
                and (not include or any(
                    [tid in include for tid in _path + [tid]]))):
            abun = counts.get(tid, 0)
        self[tid] = TaxTree(counts=abun,
                            score=scores.get(tid, NO_SCORE),
                            rank=rank,
                            acc=abun)
        if tid not in taxonomy.children:
            return self[tid].acc
        # Taxid has children (is a branch)
        for chld in taxonomy.children[tid]:
            child_acc = self[tid].allin1(
                taxonomy=taxonomy, counts=counts, scores=scores,
                ancestors=ancestors, tid=chld,
                min_taxa=min_taxa, min_rank=min_rank,
                just_min_rank=just_min_rank, include=include, exclude=exclude,
                out=out, _path=_path + [tid])
            if child_acc is None:  # No child created, continue
                continue
            rank_prune: bool = False
            # Set min_rank pruning condition using a robust algorithm
            #   suitable for multiple eukaryotic NO_RANK sublevels
            if (min_rank and (self[tid][chld].rank < min_rank or
                              rank <= min_rank or parent_rank <= min_rank)):
                rank_prune = True
            # Check conditions for pruning the leaf
            if child_acc < min_taxa or rank_prune:
                # Check conditions for saving the data of leaf in the parent
                if child_acc > 0 and (
                        not just_min_rank or (
                        rank_prune and (
                        rank == min_rank or parent_rank <= min_rank))):
                    #  and not (rank != min_rank and parent_rank > min_rank))):
                    collapsed_counts: int = (self[tid].counts
                                             + self[tid][chld].counts)
                    if collapsed_counts:  # Average the collapsed score
                        self[tid].score = swmean(
                            self[tid].counts, self[tid].score,
                            self[tid][chld].counts, self[tid][chld].score)
                        # Accumulate abundance in higher tax
                        self[tid].counts += self[tid][chld].counts
                        self[tid].acc += child_acc
                    else:  # No collapse: update acc counts erasing leaf counts
                        pass
                        # self[tid].acc -= child_acc
                pruned: TaxTree = self[tid].pop(chld)  # Prune the leaf
                assert not pruned, f'{tid}-//->{chld}->{pruned}'
                continue
            self[tid].acc += child_acc
            if out:
                populate_output(chld, self[tid])
        # If not counts, calculate score from leaves
        if not abun and self[tid].acc and self[tid]:
            self[tid].score = sum([
                self[tid][chld].score * self[tid][chld].acc
                / self[tid].acc for chld in self[tid]])
        # If ROOT, populate results
        if tid == ROOT and out:
            populate_output(ROOT, self)
        return self[tid].acc

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
                 scores: Union[Dict[TaxId, Score], SharedCounter] = None,
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
            scores: Optional I/O dict; at 1st entry should be empty.
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
                             or self[tid].rank is just_level)):
                    if abundance is not None:
                        abundance[tid] = self[tid].counts
                    if accs is not None:
                        accs[tid] = self[tid].acc
                    if scores is not None and self[tid].score != NO_SCORE:
                        scores[tid] = self[tid].score
                    if ranks is not None:
                        ranks[tid] = self[tid].rank
                if self[tid]:
                    self[tid].get_taxa(abundance, accs, scores, ranks,
                                       mindepth, maxdepth,
                                       include, exclude,
                                       just_level, in_branch)

    def shape(self) -> None:
        """
        Recursively populate accumulated counts and score.

        From bottom to top, accumulate counts in higher taxonomical
        levels, so populate self.acc of the tree. Also calculate
        score for levels that have no reads directly assigned
        (unassigned = 0). It eliminates leaves with no accumulated
        counts. With all, it shapes the tree to the most useful form.

        """
        self.acc = self.counts  # Initialize accumulated with unassigned counts
        for tid in list(self):  # Loop if this node has subtrees
            self[tid].shape()  # Accumulate for each branch/leaf
            if not self[tid].acc:
                self.pop(tid)  # Prune empty leaf
            else:
                self.acc += self[tid].acc  # Acc lower tax acc in this node
        if not self.counts:
            # If not unassigned (no reads directly assigned to the level),
            #  calculate score from leaves, trying different approaches.
            if self.acc:  # Leafs (at least 1) have accum.
                self.score = sum([self[tid].score * self[tid].acc
                                  / self.acc for tid in self])
            else:  # No leaf with unassigned counts nor accumulated
                # Currently, do nothing, but another choice is
                #   just get the averaged score by number of leaves
                # self.score = sum([self[tid].score
                #                   for tid in list(self)])/len(self)
                pass

    def subtract(self) -> None:
        """
        Recursively subtract counts of lower levels from higher ones.

        From bottom to top, subtract counts from higher taxonomical
        levels, under the assumption than they were accumulated before.

        """
        for tid in list(self):  # Loop if this node has subtrees
            self.counts -= self[tid].acc  # Subtract lower taxa acc
            self[tid].subtract()  # Subtract for each branch/leaf
            if self.counts < 0:
                self.counts = 0  # Avoid negative counts

    def prune(self,
              min_taxa: int = 1,
              min_rank: Rank = None,
              collapse: bool = True,
              debug: bool = False) -> bool:
        """
        Recursively prune/collapse low abundant taxa of the TaxTree.

        Args:
            min_taxa: minimum taxa to avoid pruning/collapsing
                one level to the parent one.
            min_rank: if any, minimum Rank allowed in the TaxTree.
            collapse: selects if a lower level should be accumulated in
                the higher one before pruning a node (do so by default).
            debug: increase output verbosity (just for debug)

        Returns: True if this node is a leaf

        """
        for tid in list(self):  # Loop if this node has subtrees
            if (self[tid]  # If the subtree has branches,
                    and self[tid].prune(min_taxa, min_rank, collapse, debug)):
                # The pruned subtree keeps having branches (still not a leaf!)
                if debug:
                    print(f'[NOT pruning branch {tid}, '
                          f'counts={self[tid].counts}]', end='')
            else:  # self[tid] is a leaf (after pruning its branches or not)
                if (self[tid].counts < min_taxa  # Not enough counts, or check
                        or (min_rank  # if min_rank is set, then check if level
                            and (self[tid].rank < min_rank  # is lower or
                                 or self.rank <= min_rank))):  # other test
                    if collapse:
                        collapsed_counts: int = self.counts + self[tid].counts
                        if collapsed_counts:  # Average the collapsed score
                            self.score = ((self.score * self.counts
                                           + self[tid].score * self[
                                               tid].counts)
                                          / collapsed_counts)
                            # Accumulate abundance in higher tax
                            self.counts = collapsed_counts
                    else:  # No collapse: update acc counts erasing leaf counts
                        if self.acc > self[tid].counts:
                            self.acc -= self[tid].counts
                        else:
                            self.acc = 0
                    if debug and self[tid].counts:
                        print(f'[Pruning branch {tid}, '
                              f'counts={self[tid].counts}]', end='')
                    self.pop(tid)  # Prune leaf
                elif debug:
                    print(f'[NOT pruning leaf {tid}, '
                          f'counts={self[tid].counts}]', end='')
        return bool(self)  # True if this node has branches (is not a leaf)

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


class MultiTree(dict):
    """Nodes of a multiple taxonomical tree"""

    def __init__(self, *args,
                 samples: List[Sample],
                 counts: Dict[Sample, int] = None,
                 accs: Dict[Sample, int] = None,
                 rank: Rank = Rank.UNCLASSIFIED,
                 scores: Dict[Sample, Score] = None
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
            scores = {sample: NO_SCORE for sample in samples}
        self.score: List[Score] = [scores[sample] for sample in samples]

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
             abundances: Dict[Sample, Counter[TaxId]] = None,
             accs: Dict[Sample, Counter[TaxId]] = None,
             scores: Dict[Sample, Dict[TaxId, Score]] = None,
             taxid: TaxId = ROOT,
             _path: List[TaxId] = None) -> None:
        """
        Recursively build a taxonomy tree.

        Args:
            taxonomy: Taxonomy object.
            abundances: Dict of counters with taxids' abundance.
            accs: Dict of counters with taxids' accumulated abundance.
            scores: Dict of dicts with taxids' score.
            taxid: It's ROOT by default for the first/base method call
            _path: list used by the recursive algorithm to avoid loops

        Returns: None

        """
        # Create dummy variables in case they are None in the 1st call
        if not _path:
            _path = []
            if not abundances:
                abundances = {sample: col.Counter({ROOT: 1})
                              for sample in self.samples}
            if not accs:
                accs = {sample: col.Counter({ROOT: 1})
                        for sample in self.samples}
            if not scores:
                scores = {sample: {} for sample in self.samples}
        if taxid not in _path:  # Avoid loops for repeated taxid (like root)
            multi_count: Dict[Sample, int] = {
                sample: abundances[sample].get(taxid, 0)
                for sample in self.samples
            }
            multi_acc: Dict[Sample, int] = {
                sample: accs[sample].get(taxid, 0)
                for sample in self.samples
            }
            multi_score: Dict[Sample, Score] = {
                sample: scores[sample].get(taxid, NO_SCORE)
                for sample in self.samples
            }
            if any(multi_acc.values()):  # Check for any populated branch
                self[taxid] = MultiTree(samples=self.samples,
                                        counts=multi_count,
                                        accs=multi_acc,
                                        scores=multi_score,
                                        rank=taxonomy.get_rank(taxid))
                if taxid in taxonomy.children:  # taxid has children
                    for child in taxonomy.children[taxid]:
                        self[taxid].grow(taxonomy=taxonomy,
                                         abundances=abundances,
                                         accs=accs,
                                         scores=scores,
                                         taxid=child,
                                         _path=_path + [taxid])

    def toxml(self,
              taxonomy: Taxonomy,
              krona: KronaTree,
              node: Elm = None,
              ) -> None:
        """
        Recursive method to generate XML.

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
                        SCORE: {self.samples[i]: (
                            f'{self[tid].score[i]:.1f}'
                            if self[tid].score[i] != NO_SCORE else '0')
                            for i in range(num_samples)},
                        }
            )
            if self[tid]:
                self[tid].toxml(taxonomy=taxonomy,
                                krona=krona,
                                node=new_node)

    def to_items(self,
                 taxonomy: Taxonomy,
                 items: List[Tuple[TaxId, List]],
                 sample_indexes: List[int] = None
                 ) -> None:
        """
        Recursive method to populate a list (used to feed a DataFrame).

        Args:
            taxonomy: Taxonomy object.
            items: Input/Output list to be populated.
            sample_indexes: Indexes of the samples of interest (for cC)

        Returns: None

        """
        for tid in self:
            list_row: List = []
            if sample_indexes:
                for i in sample_indexes:
                    list_row.append(self[tid].counts[i])
            else:
                for i in range(len(self.samples)):
                    list_row.extend([self[tid].accs[i],
                                     self[tid].counts[i],
                                     self[tid].score[i]])
                list_row.extend([taxonomy.get_rank(tid).name.lower(),
                                 taxonomy.get_name(tid)])
            items.append((tid, list_row))
            if self[tid]:
                self[tid].to_items(taxonomy=taxonomy, items=items,
                                   sample_indexes=sample_indexes)