#!/usr/bin/env python

"""
Tests for adding a data atrribute to a UGrid object

Designed to be run with pytest
"""
import numpy as np
import pytest

from pyugrid.ugrid import UGrid, DataSet
from pyugrid.test_examples import two_triangles


def test_add_all_data():
    '''
    you should not be able add a data dict directly
    '''
    grid = two_triangles()

    assert grid.data == {}

    with pytest.raises(AttributeError):
        grid.data = {'depth': DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])}


def test_add_node_data():

    grid = two_triangles()

    # create a dataset object for the depths:
    depths = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    depths.attributes['units'] = 'm'
    depths.attributes["standard_name"] = "sea_floor_depth"
    depths.attributes["positive"] = "down"

    grid.add_data(depths)

    assert grid.data['depth'].name == 'depth'
    assert grid.data['depth'].attributes['units'] == 'm'
    assert np.array_equal( grid.data['depth'].data, [1.0, 2.0, 3.0, 4.0] )

def test_add_node_data_wrong():
    """too short an array"""

    grid = two_triangles()

    # create a dataset object for the depths:
    depths = DataSet('depth', location='node', data=[1.0, 2.0, 3.0])

    with pytest.raises(ValueError):
        grid.add_data(depths)


def test_add_face_data():

    grid = two_triangles()

    # create a dataset object for velocity:
    u_vel = DataSet('u', location='face', data=[1.0, 2.0])
    u_vel.attributes['units'] = 'm/s'
    u_vel.attributes["standard_name"] = "eastward_sea_water_velocity"

    grid.add_data(u_vel)

    assert grid.data['u'].name == 'u'
    assert grid.data['u'].attributes['units'] == 'm/s'
    assert np.array_equal( grid.data['u'].data, [1.0, 2.0] )

def test_add_face_data_wrong():
    """too short an array"""

    grid = two_triangles()

    # create a dataset object for velocity:
    u_vel = DataSet('u', location='face', data=[1.0])

    with pytest.raises(ValueError):
        grid.add_data(u_vel)


def test_add_edge_data():

    grid = two_triangles()

    # create a dataset object for velocity:
    bnds = DataSet('bounds', location='edge', data=[0, 1, 0, 0, 1])
    bnds.attributes["standard_name"] = "boundary type"

    grid.add_data(bnds)

    assert grid.data['bounds'].name == 'bounds'
    assert np.array_equal( grid.data['bounds'].data, [0, 1, 0, 0, 1] )

def test_add_edge_data_wrong():
    """too long an array"""

    grid = two_triangles()

    # create a dataset object for velocity:
    # a miss-matched set
    bnds = DataSet('bounds', location='edge', data=[0, 1, 0, 0, 1, 3, 3])

    with pytest.raises(ValueError):
        grid.add_data(bnds)

def test_add_boundary_data():

    grid = two_triangles()

    print grid.boundaries

    # add the boundary definitions:
    grid.boundaries = [(0,1),
                       (0,2),
                       (1,3),
                       (2,3),
                      ]
    # create a dataset object for boundary conditions:
    bnds = DataSet('bounds', location='boundary', data=[0, 1, 0, 0, 1])
    bnds.attributes["long_name"] = "model boundary conditions"

    # wrong size for data
    with pytest.raises(ValueError):
        grid.add_data(bnds)
    # correct data
    bnds.data = [0, 1, 0, 0]
    grid.add_data(bnds)

    assert grid.data['bounds'].name == 'bounds'
    assert np.array_equal( grid.data['bounds'].data, [0, 1, 0, 0] )






