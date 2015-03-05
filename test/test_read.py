#!/usr/bin/env python

"""
Tests for testing a ugrid file read.

We really need a LOT more sample data files....

"""

import os
import contextlib
import pytest

from .utilities import chdir

import numpy as np
import netCDF4

from pyugrid import ugrid
from pyugrid import read_netcdf

UGrid = ugrid.UGrid

files = os.path.join(os.path.split(__file__)[0], 'files')
file11 = 'ElevenPoints_UGRIDv0.9.nc'

def test_simple_read():
    """
    Can it be read at all?
    """

    with chdir(files):
        ug = UGrid.from_ncfile(file11)

    assert True


def test_get_mesh_names():
    """
    Check that it can find the mesh variable.

    NOTE: this really should check for more than one mesh..
    """

    with chdir(files):
        nc = netCDF4.Dataset(file11)
    names = read_netcdf.find_mesh_names(nc)

    assert names == [u'Mesh2']


def test_mesh_not_there():
    """
    Test raising Value error with incorrect mesh name.
    """

    with pytest.raises(ValueError):
        with chdir(files):
            ug = UGrid.from_ncfile(file11, mesh_name='garbage')

def test_load_grid_from_nc():
    """
    Test reading a fairly full example file.
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11)

    assert grid.mesh_name == 'Mesh2'

    assert grid.nodes.shape == (11, 2)
    assert grid.faces.shape == (13, 3)
    assert grid.face_face_connectivity.shape == (13, 3)
    assert grid.boundaries.shape == (9, 2)

    assert grid.edges is None    # no edges in this data


def test_read_nodes():
    """
    Do we get the right nodes array?
    """

    with chdir(files):
        ug = UGrid.from_ncfile(file11)

    assert ug.nodes.shape == (11, 2)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0, :],     (-62.242, 12.774999))
    assert np.array_equal(ug.nodes[-1, :],    (-34.911235, 29.29379))


def test_read_edges():
    """
    Do we get the right edge array?
    """

    with chdir(files):
        ug = UGrid.from_ncfile(file11)

    assert ug.edges is None


def test_read_faces():
    """
    Do we get the right faces array?
    """

    with chdir(files):
        ug = UGrid.from_ncfile(file11)

    assert ug.faces.shape == (13, 3)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.faces[0, :],     (2, 3, 10))
    assert np.array_equal(ug.faces[-1, :],    (10, 5, 6))


def test_read_face_face():
    """
    Do we get the right face_face_connectivity array?
    """

    with chdir(files):
        ug = UGrid.from_ncfile(file11)

    assert ug.face_face_connectivity.shape == (13, 3)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.face_face_connectivity[0, :],     (11, 5, -1))
    assert np.array_equal(ug.face_face_connectivity[-1, :],    (-1, 5, 11))


def test_read_boundaries():
    """
    Do we get the right boundaries array?
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11)

    assert grid.boundaries.shape == (9, 2)

    # Not ideal to pull specific values out, but how else to test?
    # Note: file is 1-indexed, so these values are adjusted.
    assert np.array_equal(grid.boundaries, [[0, 1],
                                            [1, 2],
                                            [2, 3],
                                            [3, 4],
                                            [4, 0],
                                            [5, 6],
                                            [6, 7],
                                            [7, 8],
                                            [8, 5],
                                            ])


def test_read_face_coordinates():
    """
    Do we get the right face_coordinates array?
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11)

    assert grid.face_coordinates.shape == (13, 2)

    # Not ideal to pull specific values out, but how else to test?
    assert np.array_equal(grid.face_coordinates[0],
                          (-37.1904106666667, 30.57093))
    assert np.array_equal(grid.face_coordinates[-1],
                          (-38.684412, 27.7132626666667))


def test_read_edge_coordinates():
    """
    Do we get the right edge_coordinates array?
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11)

    # not in this sample file
    assert grid.edge_coordinates is None


def test_read_boundary_coordinates():
    """
    Do we get the right boundary_coordinates array?
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11)

    # not in this sample file
    assert grid.boundary_coordinates is None

def test_read_longitude_no_standard_name():

    with chdir(files):
        ug = UGrid.from_ncfile('no_stand_name_long.nc')

    assert ug.nodes.shape == (11, 2)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal(ug.nodes[0, :],     (-62.242, 12.774999))
    assert np.array_equal(ug.nodes[-1, :],    (-34.911235, 29.29379))


def test_read_data1():
    """
    Sample file has depths on the nodes -- the key should get read in.
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11, load_data=True)

    assert sorted(grid.data.keys()) == [u'boundary_count',
                                        u'boundary_types',
                                        u'depth']


def test_read_data2():
    """
    Sample file has depths on the nodes -- the data should get read in.
    """

    with chdir(files):
        grid = UGrid.from_ncfile(file11, load_data=True)

    # assert grid.data['depth'] is not None
    depth_data11 = [1, 1, 1, 102, 1, 1, 60, 1, 1, 97, 1]
    assert np.array_equal(grid.data['depth'].data, depth_data11)
    depth_attributes11 = {'standard_name': "sea_floor_depth_below_geoid",
                          'units': "m",
                          'positive': "down",
                          }
    assert grid.data['depth'].attributes == depth_attributes11

def test_read_from_nc_dataset():
    """
    minimal test, but makes sure you can read from an already open netCDF4.Dataset
    """
    with chdir(files):
        with netCDF4.Dataset(file11) as nc:
            grid = UGrid.from_nc_dataset(nc)
    
    assert grid.mesh_name == 'Mesh2'
    assert grid.nodes.shape == (11, 2)
    assert grid.faces.shape == (13, 3)


if __name__ == "__main__":
    test_simple_read()
