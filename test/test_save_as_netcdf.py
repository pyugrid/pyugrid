#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""
import numpy as np

from pyugrid.ugrid import UGrid
from pyugrid.test_examples import *

import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger('pyugrid').setLevel(logging.DEBUG)


def test_two_triangles():
    grid = two_triangles()

    grid.save_as_netcdf('2_triangles.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('2_triangles.nc')
    print ug
    
    assert ug.nodes.shape == (4,2)
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0,:], (0, 0))
    assert np.array_equal(ug.nodes[-1,:], (3, 2))
    assert np.array_equal(ug.nodes, grid.nodes)
    
    print ug.depths
    assert ug.depths.shape == (4,) 
    assert ug.depths[0] == 1 
    print ug.depth_unit
    assert ug.depth_unit == "unknown"


def test_21_triangles():
    grid = twenty_one_triangles()

    grid.save_as_netcdf('21_triangles.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('21_triangles.nc')
    print ug
    
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes, grid.nodes)
    
    print ug.depths
    assert ug.depths[0] == 1 
    print ug.depth_unit
    assert ug.depth_unit == "unknown"

if __name__ == "__main__":
    test_two_triangles()
    test_21_triangles()
