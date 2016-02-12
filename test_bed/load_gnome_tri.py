#!usr/bin/env python

"""
Test to see how to write a custom UGRID file reader.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
import netCDF4
import pyugrid


def load_gnome_tri(filename):
    """
    load an fvcom/gnome netcdf file

    :param filename: the path to the file to load

    :returns: an UGRID object with the data from the file
    """
    nc = netCDF4.Dataset(filename)

    # Check if it's the file type we are looking for.
    if nc.getncattr('grid_type').lower() != 'triangular':
        msg = ('This does not appear to be a valid triangular grid file:'
               '\nMissing the "grid_type"="Triangular" global attribute')
        raise ValueError(msg)

    # Create an empty UGrid:
    ug = pyugrid.UGrid()
    # load the nodes
    lon = nc.variables['lon']
    lat = nc.variables['lat']

    num_nodes = lon.shape[0]
    ug.nodes = np.zeros((num_nodes, 2), dtype=lon.dtype)
    ug.nodes[:, 0] = lon[:]
    ug.nodes[:, 1] = lat[:]

    # Load the faces.
    ug.faces = nc.variables['nv'][:].T

    # Load the connectivity array.
    ug.face_face_connectivity = nc.variables['nbe'][:].T

    # Load the center points of the faces.
    ug.face_coordinates = np.zeros((len(ug.faces), 2), dtype=lon.dtype)
    ug.face_coordinates[:, 0] = nc.variables['lonc'][:]
    ug.face_coordinates[:, 1] = nc.variables['latc'][:]

    bounds = nc.variables['bnd']
    # Ignoring the second two fields. What are they?
    ug.boundaries = bounds[:, :2]

    return ug


def load_from_varnames(filename, names_mapping):
    """
    load a UGrid from a netcdf file where the roles are defined by the names
    of the variables

    :param filename: the names of the file to load (or opendap url)

    :param names_mapping: dict that maps the variable names to the UGrid
    components

    """

    nc = netCDF4.Dataset("small_trigrid_example.nc")

    # check if it's the file type we are looking for
    try:
        name, value = names_mapping['attribute_check']
        if nc.getncattr(name).lower() != value:
            msg = ('This does not appear to be a valid triangular grid file:'
                   '\nMissing the "grid_type"="Triangular" global attribute')
            raise ValueError(msg)
    except KeyError:
        pass

    # Create an empty UGrid:
    ug = pyugrid.UGrid()

    # Load the nodes -- this is required
    # nodes are usually stored in two different arrays.
    lon = nc.variables[names_mapping['nodes_lon']]
    lat = nc.variables[names_mapping['nodes_lat']]

    num_nodes = lon.shape[0]
    ug.nodes = np.zeros((num_nodes, 2), dtype=lon.dtype)
    ug.nodes[:, 0] = lon[:]
    ug.nodes[:, 1] = lat[:]

    # Load the faces.
    faces = nc.variables[names_mapping['faces']]
    # Does it need to be transposed?
    # Assume there are more than three triangles.
    if faces.shape[0] <= faces.shape[1]:
        # Fortran order -- needs to be transposed.
        faces = faces[:].T
    else:
        faces = faces[:]
    # Is it one-indexed?
    if faces.min() == 1:
        one_indexed = True
        faces -= 1
        ug.faces = faces
    else:
        one_indexed = False

    # Load the connectivity array: optional.
    if 'face_face_connectivity' in names_mapping:
        face_face_connectivity = nc.variables[names_mapping['face_face_connectivity']]  # noqa
        # Does it need to be transposed?
        # Assume there are more than three triangles.
        if face_face_connectivity.shape[0] <= face_face_connectivity.shape[1]:
            # Fortran order -- needs to be transposed.
            face_face_connectivity = face_face_connectivity[:].T
        else:
            face_face_connectivity = face_face_connectivity[:]
        if one_indexed:
            face_face_connectivity -= 1
        ug.face_face_connectivity = face_face_connectivity[:]

    # Load the center points of the faces: optional.
    if 'face_coordinates_lon' in names_mapping and 'face_coordinates_lon' in names_mapping:  # noqa
        ug.face_coordinates = np.zeros((len(ug.faces), 2), dtype=lon.dtype)
        ug.face_coordinates[:, 0] = nc.variables[names_mapping['face_coordinates_lon']][:]  # noqa
        ug.face_coordinates[:, 1] = nc.variables[names_mapping['face_coordinates_lat']][:]  # noqa

    print("checking bounds")
    if 'boundaries' in names_mapping:  # Optional.
        print("loading boundaries")
        # FIXME: this one is weird and non-conforming!
        boundaries = nc.variables[names_mapping['boundaries']][:]
        if one_indexed:
            boundaries[:, :2] -= 1
        # Ignoring the rest of the fields. They are boundary data (see below).
        ug.boundaries = boundaries[:, :2]
        # Kludge to handle the gnome tri format's putting boundary data in
        # with the boundary indexes.
        # This looks like GNOME-style boundary data.
        if boundaries.shape[1] == 4:
            # Add the boundary data as uvars.
            pass
    return ug


if __name__ == "__main__":

    names_mapping = {
        'attribute_check': ('grid_type', 'triangular'),
        'nodes_lon': 'lon',
        'nodes_lat': 'lat',
        'faces': 'nv',
        'face_face_connectivity': 'nbe',
        'face_coordinates_lon': 'lonc',
        'face_coordinates_lat': 'latc',
        'boundaries': 'bnd',
        'associated_variables': [('u', 'face'),
                                 ('v', 'face'),
                                 ('a1u', 'face'),
                                 ('a2u', 'face'),
                                 ],
    }

    ug = load_from_varnames("small_trigrid_example.nc", names_mapping)

    print(ug)
