# pyugrid code
import pyugrid.ugrid

# test stuff
from write_nc_test_files import quad_and_triangle


def test_read_flexible_mesh(quad_and_triangle):
    """
    Test if we get back a mesh from a flexible mesh
    """
    grid = pyugrid.ugrid.UGrid.from_nc_dataset(quad_and_triangle)
    assert grid.mesh_name == 'Mesh2'

def test_read_flexible_mesh_mask(quad_and_triangle):
    """
    Test if we get back a masked array from a flexible mesh (for faces and edges)
    """
    grid = pyugrid.ugrid.UGrid.from_nc_dataset(quad_and_triangle)
    assert grid.mesh_name == 'Mesh2'
    assert hasattr(grid.faces, 'mask'), "expected masked faces"


def test_read_flexible_mesh_nodes_per_face(quad_and_triangle):
    """
    Test if we the grid contains both triangles and quads
    """
    grid = pyugrid.ugrid.UGrid.from_nc_dataset(quad_and_triangle)
    n_nodes_per_face = (~grid.faces.mask).sum(axis=1)
    assert set(n_nodes_per_face) == set([3, 4]), 'expected triangles and quads'
