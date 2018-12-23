"""
Taxonomy class, currently representing the Gene Ontology.

"""
from abc import ABCMeta, abstractmethod
from typing import Set, Union, Iterable, Tuple

from recentrifuge.config import Id, Parents, Children
from recentrifuge.rank import Rank


class Ontology(metaclass=ABCMeta):
    """Ontology abstract base class"""

    # Abstract attributes
    ROOT: Id
    parents: Parents
    children: Children

    # Default attributes
    collapse: bool = True
    excluding: Union[Tuple, Set[Id]] = ()
    including: Union[Tuple, Set[Id]] = ()
    debug: bool = False

    @abstractmethod
    def get_rank(self, anid: Id) -> Rank:
        """Retrieve the rank"""
        raise NotImplementedError(f'Class {self.__class__.__name__} does not'
                                  f' implement get_rank()')

    @abstractmethod
    def get_name(self, anid: Id) -> str:
        """Retrieve the name"""
        raise NotImplementedError(f'Class {self.__class__.__name__} does not'
                                  f' implement get_name()')

    @abstractmethod
    def get_ancestors(self, leaves: Iterable
                      ) -> Tuple[Set, Set]:
        """Return the ids entered with all their ancestors"""
        raise NotImplementedError(f'Class {self.__class__.__name__} does not'
                                  f' implement get_ancestors()')
