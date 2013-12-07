#!/usr/bin/env python

"""
tests for testing a ugrid file read.

we really need some  more sample data files....

"""
import numpy as np

from pyugrid.ugrid import UGrid


def test_simple_read():
	""" can it be read at all """
	ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

	assert True

def test_read_nodes():
	""" Do we get the right nodes array? """
	ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

	assert ug.nodes.shape == (11,2)

	# not ideal to pull specific values out, but how else to test?
	assert np.array_equal( ug.nodes[0,:],	 (-62.242, 12.774999) )
	assert np.array_equal( ug.nodes[-1,:],	 (-34.911235,  29.29379) )

## no edge data in test file at this point
# def test_read_edges():
# 	""" Do we get the right edges array? """
# 	ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

# 	print ug.edges

# 	assert False

def test_read_face_node_connectivity():
	""" Do we get the right connectivity array? """
	ug = UGrid.from_ncfile('files/ElevenPoints_UGRIDv0.9.nc')

	assert ug.faces.shape == (13, 3)

	# # not ideal to pull specific values out, but how else to test?
	## note: file is 1-indexed, so these values are adjusted
	assert np.array_equal( ug.faces[0,:],	 (2, 3, 10) )
	assert np.array_equal( ug.faces[-1,:],	 (10, 5, 6) )

# def test_simple_read():
# 	ug = UGrid.from_ncfile('files/two_triangles.nc')

# 	assert False

 