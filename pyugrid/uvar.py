#!/usr/bin/env python

"""
UVar object, used to hold variables that are associated with a ugrid
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from .util import asarraylike

class UVar(object):
    """
    A class to hold a variable associated with the UGrid. Data can be on the 
    nodes, edges, etc. -- "UGrid Variable"

    It holds an array of the data, as well as the attributes associated
    with that data  -- this is mapped to a netcdf variable with 
    attributes(attributes get stored in the netcdf file)
    """

    def __init__(self, name, location='node', data=None, attributes=None):
        """
        create a UVar object
        :param name: the name of the data (depth, u_velocity, etc.)
        :type name: string

        :param location: the type of grid element the data is associated with:
                         'node', 'edge', or 'face' the data is assigned to

        :param data: the data
        :type data: 1-d numpy array or array-like object ().
                    IF you have a list or tuple, it should be converted or something compatible (list, etc.)        
        """
        self.name = name

        if location not in ['node', 'edge', 'face', 'boundary']:
            raise ValueError("location must be one of: 'node', 'edge', 'face', 'boundary'")
 
        self.location = location 

        if data is None:
            self._data = np.zeros((0,), dtype=np.float64) # could be any data type
        else:
            self._data = asarraylike(data)

        self.attributes = {} if attributes is None else attributes

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, data):
        self._data = asarraylike(data)
    @data.deleter
    def data(self):
        self._data = self._data = np.zeros((0,), dtype=np.float64)

    def __str__(self):
        print ("in __str__, data is:", self.data)
        return "UVar object: {0:s}, on the {1:s}s, and {2:d} data points\nAttributes: {3}".format(self.name, self.location, len(self.data), self.attributes)


