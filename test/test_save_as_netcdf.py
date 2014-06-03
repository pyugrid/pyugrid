#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""
import numpy as np

from pyugrid.ugrid import UGrid
from pyugrid.test_examples import two_triangles

import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger('pyugrid').setLevel(logging.DEBUG)


def test_simple_write():

    grid = two_triangles()
    grid.depths = [1,2,3,4]
    grid.depth_unit = "unknown"

    grid.save_as_netcdf('two_triangles.nc')

    ## be good to have an actual test here...
    assert True


def test_simple_read():
    """ Do we get the right nodes array? """
    ug = UGrid.from_ncfile('two_triangles.nc')
    print ug
    
    assert ug.nodes.shape == (4,2)
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal( ug.nodes[0,:],	 (0, 0) )
    assert np.array_equal( ug.nodes[-1,:],	 (3, 2) )
    
    print ug.depths
    assert ug.depths.shape == (4,) 
    assert ug.depths[0] == 1 
    print ug.depth_unit
    assert ug.depth_unit == "unknown"

if __name__ == "__main__":
    test_simple_write()
    test_simple_read()
