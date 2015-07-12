#!/usr/bin/env python

"""
Tests for the UVar object
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
import pytest

from pyugrid.ugrid import UVar


def test_init():
    # create a UVar object for the depths:
    d = UVar('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    assert d.name == 'depth'
    print( d )
    print ( d.data )
    assert np.array_equal( d.data, [1.0, 2.0, 3.0, 4.0] )
    assert d.location == 'node'

    assert d.attributes == {}

    with pytest.raises(ValueError):       
        d = UVar('depth', location='nodes')

def test_add_data():

    d = UVar('depth', location='node')

    assert d.name == 'depth'
    assert np.array_equal(d.data, [] )

    #add the data:
    d.data = [1.0, 2.0, 3.0, 4.0]

    assert np.array_equal(d.data, [1.0, 2.0, 3.0, 4.0] )
    # duck type check of ndarray
    d.data *= 2
    assert np.array_equal(d.data, [2.0, 4.0, 6.0, 8.0] )

def test_delete_data():
    # create a UVar object for the depths:
    d = UVar('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    del d.data

    assert np.array_equal(d.data, [] )

def test_str():
    d = UVar('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])
    print(str(d))
    assert str(d) == "UVar object: depth, on the nodes, and 4 data points\nAttributes: {}"

def add_attributes():
    d = UVar('depth', location='node', data=[1.0, 2.0, 3.0, 4.0])

    d.attributes = {"standard_name" : "sea_floor_depth_below_geoid",
                    "units" : "m",
                    "positive" : "down",
                    }

    assert d.attributes['units'] == 'm'
    assert d.attributes['posative'] == 'down'



 