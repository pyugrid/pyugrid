#!/usr/bin/env python

"""
testing of code to find nodes:
  currently only nearest neighbor

"""

import numpy as np

from pyugrid.test_examples import two_triangles, twenty_one_triangles

def test_locate_node():
    """
    test finding a single node
    """

    ugrid = twenty_one_triangles()

    assert ugrid.locate_nodes( (4.58, 5.08) ) == 6

def test_locate_nodes():
    """
    test finding multiple nodes at once
    """

    ugrid = twenty_one_triangles()

    assert np.array_equal( ugrid.locate_nodes( ((4.58, 5.08),
                                                (4.81, 0.89),
                                                (6.43, 12.9),
                                                (8.74, 6.86),
                                                (5.12, 7.31),
                                                )),
                            (6, 0, 17, 10, 8)
                          )

def test_locate_exact():
    """
    The nearest neighbor of the exact node locations has better be the nodes!
    """
    ugrid = twenty_one_triangles()

    assert np.array_equal( ugrid.locate_nodes( ugrid.nodes ),
                        range(len(ugrid.nodes))
                       )    


if __name__ == "__main__":
    test_locate_nodes()


