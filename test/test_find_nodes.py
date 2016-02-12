#!/usr/bin/env python

"""
testing of code to find nodes:
  currently only nearest neighbor

"""

import numpy as np
import pytest

from pyugrid.test_examples import twenty_one_triangles


@pytest.mark.slow
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
    The nearest neighbor of the exact node locations had better be
    the nodes!
    """
    ugrid = twenty_one_triangles()

    assert np.array_equal( ugrid.locate_nodes( ugrid.nodes ),
                        range(len(ugrid.nodes))
                       )    

def test_locate_middle():
    """
    see what happens the point is equidistant to two nodes
    """
    ugrid = twenty_one_triangles()

    # (3,5) is equidistant bewteen nodes 2, 6 and 7
    # 2 is returned, but might be arbitrary
    # assert ugrid.locate_nodes( (3, 5) ) == 2

    # perturb the point a bit, and nearest changes to:
    assert ugrid.locate_nodes( (3.0000000001, 5)  ) == 6
    assert ugrid.locate_nodes( (3, 5.00000000001) ) == 7
    assert ugrid.locate_nodes( (3, 4.99999999999) ) == 2


if __name__ == "__main__":
    test_locate_nodes()


