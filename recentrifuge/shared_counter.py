"""
Contains class useful with operations with shared counters of taxa
"""

import collections as col
# from typing import Counter


class SharedCounter(col.Counter):
    """Extends collection.Counter with useful ops. for shared taxa."""
    # Disable warning of 'fromkeys' not implemented, as it is yet not
    #  implemented in Counter class.
    # pylint: disable=abstract-method

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

    def __mul__(self, other):
        """c * d multiply each element of c by the element in d, if exists."""
        if isinstance(self, col.Counter) and isinstance(other, col.Counter):
            result: SharedCounter = SharedCounter()
            for item in self:
                if item in other:
                    result[item] = self[item] * other[item]
            return result
        return NotImplemented

    def __imul__(self, other):
        """c *= d multiply each element of c by the element in d, if exists."""
        self = SharedCounter.__mul__(self, other)
        return self

    def __truediv__(self, other):
        """c / d divide each element of c by the element in d, if exists."""
        if isinstance(self, col.Counter) and isinstance(other, col.Counter):
            result: SharedCounter = SharedCounter()
            for item in self:
                if item in other:
                    result[item] = self[item] / other[item]  # type: ignore
            return result
        return NotImplemented

    def __itruediv__(self, other):
        """c /= d divide each element of c by the element in d, if exists."""
        self = SharedCounter.__truediv__(self, other)
        return self

    def __floordiv__(self, other):
        """c // i floor divide each element of c by integer i."""
        if isinstance(self, col.Counter) and isinstance(other, int):
            result: SharedCounter = SharedCounter()
            for item in self:
                result[item] = self[item] // other
            return result
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

    def __pos__(self):
        """+d just remove non-positive counts."""
        return SharedCounter(col.Counter.__pos__(self))
