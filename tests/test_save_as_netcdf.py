#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""

from pyugrid.ugrid import UGrid


# a really simple grid:

def two_triangles():
    """
    returns about the simplest triangle grid possible
    """
    nodes = [(0.0, 0.0),
             (2.0,0.0),
             (1.0,2.0),
             (3.0,2.0)]

    faces = [(0, 1, 2),
             (1, 3, 2),]

    edges = [(0,1),
             (1,3),
             (3,2),
             (2,0)]
    return UGrid(nodes, faces, edges)


def test_simple_write():

    grid = two_triangles()

    grid.save_as_netcdf('temp/two_triangle.nc')

    assert False

