#!/usr/bin/env python

"""
tests for saving a UGrid in netcdf format

designed to be run with pytest
"""

from pyugrid.ugrid import UGrid
from pyugrid.test_examples import two_triangles


def test_simple_write():

    grid = two_triangles()

    grid.save_as_netcdf('two_triangles.nc')

    ## be good to have an actual test here...
    assert True

 