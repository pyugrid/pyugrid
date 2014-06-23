#!/usr/bin/env python

"""
testing of various utilities to manipulate the grid

"""

import numpy as np

from pyugrid.test_examples import two_triangles, twenty_one_triangles

def test_build_face_face_connectivity():
    ugrid = two_triangles()

    ugrid.build_face_face_connectivity()

    face_face = ugrid.face_face_connectivity

    #zeroth triangle should have two empty faces
    assert sum(face_face[0] == -1) == 2
    #zeroth triangle should have one neighbor: 1
    assert sum(face_face[0] == 1) == 1
    #first triangle should have two empty faces
    assert sum(face_face[1] == -1) == 2
    #first triangle should have one neighbor: 0
    assert sum(face_face[1] == 1) == 0

def test_build_face_face_connectivity2():
    """
    test with a slightly larger mesh
    """
    ugrid = twenty_one_triangles()

    ugrid.build_face_face_connectivity()

    face_face = ugrid.face_face_connectivity

    assert np.sort(face_face[0]).tolist() ==  [-1,  2,  3]
    assert np.sort(face_face[8]).tolist() ==  [-1,  6,  9]
    assert np.sort(face_face[15]).tolist() == [13, 14, 16]
    assert np.sort(face_face[20]).tolist() == [-1, -1, 19]
    assert np.sort(face_face[9]).tolist() ==  [ 7,  8, 10]

def test_build_edges():
    ugrid = two_triangles()
    ugrid.build_edges()
    edges = ugrid.edges

    edges.sort(axis=0)
    print edges
    assert np.array_equal(edges, [[0, 1],
                                  [0, 2],
                                  [1, 2],
                                  [1, 3],
                                  [2, 3]])

def test_build_face_coordinates():
    grid = two_triangles()
    grid.build_face_coordinates()
    coords = grid.face_coordinates

    assert coords.shape == (2,2)
    assert np.allclose(coords, [ (1.1, 0.76666667),
                                 (2.1, 1.43333333)]) 

def test_build_edge_coordinates():
    grid = two_triangles()
    grid.build_edge_coordinates()
    coords = grid.edge_coordinates
    nodes = grid.nodes
    print coords

    assert coords.shape == (5,2)
    assert np.allclose(coords, [[ 1.1,  0.1],
                                [ 2.6,  1.1],
                                [ 2.1,  2.1],
                                [ 0.6,  1.1],
                                [ 1.6,  1.1]]) 

def test_build_boundary_coordinates():
    grid = two_triangles()
    # add some boundaries
    grid.boundaries = [(0,1),
                       (0,2),
                       (2,3),
                       (1,3)]
    grid.build_boundary_coordinates()
    coords = grid.boundary_coordinates
    nodes = grid.nodes
    print coords

    assert coords.shape == (4,2)
    assert np.allclose(coords, [[ 1.1, 0.1],
                                [ 0.6, 1.1],
                                [ 2.1, 2.1],
                                [ 2.6, 1.1]]) 





