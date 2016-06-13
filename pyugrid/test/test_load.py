#!/usr/bin/env python

"""
tests for pyugrid.load()

not much here, but maybe there will be more logic later
"""
import os
import netCDF4
import pyugrid
from .utilities import chdir

files = os.path.join(os.path.split(__file__)[0], 'files')


def test_load_filename():
    """
    load from a filename
    """
    with chdir(files):
        ug = pyugrid.load('ElevenPoints_UGRIDv0.9.nc')
    assert isinstance(ug, pyugrid.UGrid)


def test_load_dataset():
    """
    load from an open Dataset
    """
    with chdir(files):
        ds = netCDF4.Dataset('ElevenPoints_UGRIDv0.9.nc')
        ug = pyugrid.load(ds)
    assert isinstance(ug, pyugrid.UGrid)
