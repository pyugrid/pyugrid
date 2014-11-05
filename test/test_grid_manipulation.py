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

    # order matters!
    print ugrid.faces
    assert np.array_equal( face_face[0], [-1, 1, -1] )
    assert np.array_equal( face_face[1], [-1, -1, 0] )

def test_build_face_face_connectivity2():
    """
    test with a slightly larger mesh
    """
    ugrid = twenty_one_triangles()

    ugrid.build_face_face_connectivity()

    face_face = ugrid.face_face_connectivity

    assert face_face[0].tolist()  ==  [-1,  3,  2]
    assert face_face[8].tolist()  ==  [-1,  9,  6]
    assert face_face[15].tolist() ==  [14, 16, 13]
    assert face_face[20].tolist() ==  [19, -1, -1]
    assert face_face[9].tolist()  ==  [ 8,  10, 7]

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


def test_build_boundaries():
    ugrid = two_triangles()

    ugrid.build_face_face_connectivity()
    ugrid.build_boundaries()

    boundaries = ugrid.boundaries.tolist()
    boundaries.sort() # it doesn' matter what order they are in

    assert boundaries == [[0, 1], [1, 3], [2, 0], [3, 2]]
    
def test_build_boundaries2():
    """
    same as above, but with larger set
    """
    ugrid = twenty_one_triangles()

    ugrid.build_face_face_connectivity()
    ugrid.build_boundaries()

    boundaries = ugrid.boundaries.tolist()
    boundaries.sort() # it doesn't matter what order they are in

    assert boundaries == [[0, 1], [1, 5], [2, 0], [3, 6], [4, 3], [5, 11], [6, 9], [7, 2], [9, 10], [10, 4], [11, 14], [12, 7], [13, 12], [14, 16], [15, 13], [16, 18], [17, 15], [18, 19], [19, 17]]
    

