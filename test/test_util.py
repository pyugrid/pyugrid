"""
Tests for the utilities module.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from pyugrid import util


class DummyArrayLike():
    """
    Class that will look like an array to this function, even
    though it won't work!

    Just for tests.

    This will need to be updated when the function is changed.

    """
    def dtype(self):
        pass

    def shape(self):
        pass


def test_asarraylike_list():
    """
    Passing in a list should return a np.ndarray.

    """
    lst = [1, 2, 3, 4]
    result = util.asarraylike(lst)
    assert isinstance(result, np.ndarray)
    assert np.array_equal(result, lst)


def test_asarraylike_array():
    """
    Passing in a list should return a np.ndarray.

    """
    arr = np.array([1, 2, 3, 4])
    result = util.asarraylike(arr)

    assert result is arr


def test_as_test_asarraylike_dummy():
    dum = DummyArrayLike()
    result = util.asarraylike(dum)
    assert result is dum
