#!/usr/bin/env python

"""
tests for finding a given DataSet by standard name

This could use some more testing and cleaning up -- no need for all
this reading and writting -- that is tested elsewhere

designed to be run with pytest
"""
import numpy as np

from pyugrid.ugrid import UGrid, DataSet
from pyugrid.test_examples import *

import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger('pyugrid').setLevel(logging.DEBUG)


def two_triangles_with_depths():
    grid = two_triangles()

    depths = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    depths.attributes['units'] = 'unknown'
    depths.attributes['standard_name'] = "sea_floor_depth_below_geoid"
    depths.attributes['positive'] = "down"
    grid.add_data(depths)

    return grid

def find_depths(grid):
    found = grid.find_data_sets('sea_floor_depth_below_geoid')
    if found:
        return found.pop()
    return None

def test_no_std_name():
    """
    tests to make sure it doesn' crash if a DataSet does not have a standard_name
    """
    grid = two_triangles_with_depths()

    junk = DataSet('junk', location='node', data=[1.0, 2.0, 3.0, 4.0])
    junk.attributes['units'] = 'unknown'
    grid.add_data(junk)
    
    depths = find_depths(grid)

    assert depths.name == 'depth'

    
def test_two_triangles():
    grid = two_triangles_with_depths()

    grid.save_as_netcdf('2_triangles.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('2_triangles.nc', load_data=True)
    
    assert ug.nodes.shape == (4,2)
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0,:], (0.1, 0.1))
    assert np.array_equal(ug.nodes[-1,:], (3.1, 2.1))
    assert np.array_equal(ug.nodes, grid.nodes)
    
    depths = find_depths(ug)
    assert depths.data.shape == (4,) 
    assert depths.data[0] == 1
    assert depths.attributes['units'] == "unknown"

def twenty_one_triangles_with_depths():
    """
    returns a basic triangle grid with 21 triangles, a hole and a "tail"
    """
    grid = twenty_one_triangles()

    depths = DataSet('depth', location='node', data=range(1, 21))
    depths.attributes['units'] = 'unknown'
    depths.attributes['standard_name'] = "sea_floor_depth_below_geoid"
    depths.attributes['positive'] = "down"
    grid.add_data(depths)

    return grid

def test_21_triangles():
    grid = twenty_one_triangles_with_depths()

    grid.save_as_netcdf('21_triangles.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('21_triangles.nc', load_data=True)
    
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes, grid.nodes)
    
    depths = find_depths(ug)
    assert depths.data.shape == (20,) 
    assert depths.data[0] == 1
    assert depths.attributes['units'] == "unknown"

def test_two_triangles_without_faces():
    grid = two_triangles_with_depths()
    grid.faces = None

    grid.save_as_netcdf('2_triangles_without_faces.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('2_triangles_without_faces.nc', load_data=True)
    
    assert ug.nodes.shape == (4,2)
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0,:], (0.1, 0.1))
    assert np.array_equal(ug.nodes[-1,:], (3.1, 2.1))
    assert np.array_equal(ug.nodes, grid.nodes)
    
    assert ug.faces is None
    
    assert ug.edges.shape == grid.edges.shape
    assert np.array_equal(ug.edges[0,:], (0, 1))
    assert np.array_equal(ug.edges[3,:], (2, 0))
    
    depths = find_depths(ug)
    assert depths.data.shape == (4,) 
    assert depths.data[0] == 1
    assert depths.attributes['units'] == "unknown"

def test_two_triangles_without_edges():
    grid = two_triangles_with_depths()
    grid.edges = None

    grid.save_as_netcdf('2_triangles_without_edges.nc')

    # read it back in and check it out
    ug = UGrid.from_ncfile('2_triangles_without_edges.nc', load_data=True)
    
    assert ug.nodes.shape == (4,2)
    assert ug.nodes.shape == grid.nodes.shape
    
    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0,:], (0.1, 0.1))
    assert np.array_equal(ug.nodes[-1,:], (3.1, 2.1))
    assert np.array_equal(ug.nodes, grid.nodes)
    
    assert ug.faces.shape == grid.faces.shape
    
    assert ug.edges is None
    
    depths = find_depths(ug)
    assert depths.data.shape == (4,) 
    assert depths.data[0] == 1
    assert depths.attributes['units'] == "unknown"


if __name__ == "__main__":
    test_two_triangles()
    test_21_triangles()
    test_two_triangles_without_faces()
    test_two_triangles_without_edges()
