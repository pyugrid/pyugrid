import numpy as np

from pyugrid.test_examples import twenty_one_triangles

def test_single_simple():
    ugrid = twenty_one_triangles()
    face = ugrid.locate_faces((4,6.5),True)
    assert face == 6

def test_multi_simple():
    ugrid = twenty_one_triangles()
    face = ugrid.locate_faces(np.array(((4,6.5),(7,2))), True)
    assert (face == np.array((6,0))).all()

def test_single_celltree():
    ugrid = twenty_one_triangles()
    face = ugrid.locate_faces((4,6.5))
    assert face == 6

def test_multi_celltree():
    ugrid = twenty_one_triangles()
    face = ugrid.locate_faces(np.array(((4,6.5),(7,2))))
    assert (face == np.array((6,0))).all()

def test_oob():
    ugrid = twenty_one_triangles()
    face = ugrid.locate_faces((0,0),True)
    assert face == -1
    face = 0
    face = ugrid.locate_faces(np.array(((0,0),)))
    assert face == np.array((-1))