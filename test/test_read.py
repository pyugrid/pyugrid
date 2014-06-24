#!/usr/bin/env python

"""
Tests for testing a ugrid file read.

We really need a LOT more sample data files....

"""
import pytest

import numpy as np
import netCDF4

from pyugrid import ugrid
from pyugrid import read_netcdf

UGrid = ugrid.UGrid

def test_simple_read():
    """ can it be read at all """
    ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert True

def test_get_mesh_names():
    """
    check that it can find the mesh variable
    
    NOTE: this really should check for more than one mesh..
    """
    nc = netCDF4.Dataset('files/ElevenPoints_UGRIDv0.9.nc')
    names = read_netcdf.find_mesh_names( nc )

    assert names == [u'Mesh2']

def test_mesh_not_there():
    with pytest.raises(ValueError):
        ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc', mesh_name='garbage')

def test_load_grid_from_nc():
    """
    test reading a fairly full example file
    """
    grid = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert grid.mesh_name == 'Mesh2'

    assert grid.nodes.shape == (11, 2)
    assert grid.faces.shape == (13, 3)
    assert grid.face_face_connectivity.shape == (13,3) 
    assert grid.boundaries.shape == (9,2)

    assert grid.edges is None # no edges in this data

def test_read_nodes():
	""" Do we get the right nodes array? """
	ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

	assert ug.nodes.shape == (11,2)

	# not ideal to pull specific values out, but how else to test?
	assert np.array_equal( ug.nodes[0,:],	 (-62.242, 12.774999) )
	assert np.array_equal( ug.nodes[-1,:],	 (-34.911235,  29.29379) )

def test_read_edges():
    """ sample file has no edges """
    ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert ug.edges is None

def test_read_faces():
    """ Do we get the right faces array? """
    ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert ug.faces.shape == (13, 3)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal( ug.faces[0,:],    ( 2, 3, 10) )
    assert np.array_equal( ug.faces[-1,:],   (10, 5,  6) )

def test_read_face_face():
    """ Do we get the right face_face_connectivity array? """
    ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert ug.face_face_connectivity.shape == (13, 3)

    # not ideal to pull specific values out, but how else to test?
    assert np.array_equal( ug.face_face_connectivity[0,:],    ( 11, 5, -1) )
    assert np.array_equal( ug.face_face_connectivity[-1,:],   (-1, 5,  11) )



def test_read_boundaries():
    """ Do we get the right boundaries array? """
    grid = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert grid.boundaries.shape == (9, 2)

    # # not ideal to pull specific values out, but how else to test?
    ## note: file is 1-indexed, so these values are adjusted
    assert np.array_equal( grid.boundaries,[[0, 1],
                                            [1, 2],
                                            [2, 3],
                                            [3, 4],
                                            [4, 0],
                                            [5, 6],
                                            [6, 7],
                                            [7, 8],
                                            [8, 5]]
                                            )

def test_read_face_coordinates():
    """ Do we get the right face_coordinates array? """
    grid = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    assert grid.face_coordinates.shape == (13, 2)

    # # not ideal to pull specific values out, but how else to test?
    assert np.array_equal( grid.face_coordinates[0], (-37.1904106666667, 30.57093) )
    assert np.array_equal( grid.face_coordinates[-1], (-38.684412, 27.7132626666667) )

def test_read_edge_coordinates():
    grid = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    # not in this sample file
    assert grid.edge_coordinates is None

def test_read_boundary_coordinates():
    """ Do we get the right boundary_coordinates array? """
    grid = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

    # not in this sample file
    assert grid.boundary_coordinates is None


if __name__ == "__main__":
    test_simple_read()

 