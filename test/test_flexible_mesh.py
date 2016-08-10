# builtins
import os
import logging

# external
import pytest

# pyugrid code
import pyugrid.ugrid

# test stuff
from write_nc_test_files import quad_and_triangle


def test_read_flexible_mesh_nodes(quad_and_triangle):
    """
    Test if we get back nodes from a flexible mesh

    """
    grid = pyugrid.ugrid.UGrid.from_nc_dataset(quad_and_triangle)
    assert grid.mesh_name == 'Mesh2'
    assert hasattr(grid.faces, 'mask'), "we should get back a masked array"
