#!/usr/bin/env python

"""
Simple test code for getting names columns and np.ndarray.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np


class LonLatArray(np.ndarray):

    def __new__(cls, data):
        """
        Custom new,  so that we can call the "array" functionality.

        """
        arr = np.asarray(data, dtype=np.float64)
        if arr.shape[1] != 2:
            raise ValueError("data must be Nx2 array")

        return arr.view(cls)

    @property
    def lon(self):
        return self[:, 0]

    @property
    def lat(self):
        return self[:, 1]


if __name__ == "__main__":

    # Create one, just like a regular array:
    nodes = LonLatArray([(-123.02, 45.22),
                         (-123.01, 45.23),
                         (-123.03, 45.24),
                         (-123.02, 45.03),
                         (-123.1,  45.00),
                         (-123.08, 44.89),
                         ])

    # Access the whole thing or the individual lon and lat arrays.
    print(nodes)
    print(nodes.lat)
    print(nodes.lon)

    # Do math on it:
    nodes += (-1.0, 2.0)
    print(nodes)

    # or on part:
    nodes.lon[:] *= 1.1
    print(nodes)
