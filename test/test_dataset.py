#!/usr/bin/env python

"""
Tests for the DataSet object
"""

import numpy as np
import pytest

from pyugrid.ugrid import DataSet


def test_init():
    # create a dataset object for the depths:
    d = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    assert d.name == 'depth'
    assert np.array_equal( d.data, [1.0, 2.0, 3.0, 4.0] )
    assert d.location == 'node'

    assert d.attributes == {}

    with pytest.raises(ValueError):       
        d = DataSet('depth', location='nodes')

def test_add_data():

    d = DataSet('depth', location='node')

    assert d.name == 'depth'
    assert np.array_equal(d.data, [] )

    #add the data:
    d.data = [1.0, 2.0, 3.0, 4.0]

    assert np.array_equal(d.data, [1.0, 2.0, 3.0, 4.0] )
    # duck type check of ndarray
    d.data *= 2
    assert np.array_equal(d.data, [2.0, 4.0, 6.0, 8.0] )

def test_delete_data():
    # create a dataset object for the depths:
    d = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    del d.data

    assert np.array_equal(d.data, [] )

def test_str():
    d = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    print str(d)
    assert str(d) == "DataSet object: depth, on the nodes, and 4 data points\nAttributes: {}"

def add_attributes():
    d = DataSet('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    d.attributes = {"standard_name" : "sea_floor_depth_below_geoid",
                    "units" : "m",
                    "positive" : "down",
                    }

    assert d.attributes['units'] == 'm'
    assert d.attributes['posative'] == 'down'



 