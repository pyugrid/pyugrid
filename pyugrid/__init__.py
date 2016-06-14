"""
__init__.py for pyugrid package

This brings in the names we want in the package.

"""

from __future__ import (absolute_import, division, print_function)

from .ugrid import UGrid
from .uvar import UVar
from .uvar import UMVar
from . import grid_io

__version__ = '0.1.8'

__all__ = ['UGrid', 'UVar', 'UMVar', 'grid_io', 'load']


def load(filename, load_data=False):
    """
    load a UGRid object from one of:

    - open netCDF4 Dataset object
    - filename
    - OpenDAP url

    :param filename: name of file, or URL or open netCDF4 Dataset object

    :returns: The loaded UGrid object
    """
    import netCDF4

    if isinstance(filename, netCDF4.Dataset):
        return UGrid.from_nc_dataset(filename, load_data=load_data)
    else:
        return UGrid.from_ncfile(filename, load_data=load_data)
