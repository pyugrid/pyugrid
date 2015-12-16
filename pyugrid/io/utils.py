#!/usr/bin/env python

"""
utilities to help with io
"""

from __future__ import (absolute_import, division, print_function)

def load_from_varnames(filename, names_mapping):
    """
    load a UGrid from a netcdf file where the roles are defined by the names of the variables

    :param filename: the names of the file to load (or opendap url)

    :param names_mapping: dict that maps the variable names to the UGRid components

    """

    nc = netCDF4.Dataset("small_trigrid_example.nc")

    # check if it's the file type we are looking for
    try:
        name, value = names_mapping['attribute_check']
        if nc.getncattr(name).lower() != value:
            raise ValueError('This does not appear to be a valid triangular grid file:\nIt does not have the "grid_type"="Triangular" global attribute')
    except KeyError:
        pass

    # create an empty UGrid:
    ug = pyugrid.UGrid()

    # load the nodes -- this is required
    # nodes are usually stored in two different arrays
    lon = nc.variables[names_mapping['nodes_lon']]
    lat = nc.variables[names_mapping['nodes_lat']]

    num_nodes = lon.shape[0]
    ug.nodes = np.zeros((num_nodes, 2), dtype=lon.dtype)
    ug.nodes[:,0] = lon[:]
    ug.nodes[:,1] = lat[:]

    # load the faces
    faces = nc.variables[names_mapping['faces']]
    # does it need to be transposed?
    if faces.shape[0] <= faces.shape[1]:  # assume there are more than three triangles...
        # fortran order -- needs to be transposed
        faces = faces[:].T
    else:
        faces = faces[:]
    # is it one-indexed?
    if faces.min() == 1:
        one_indexed = True
        faces -= 1
        ug.faces = faces
    else:
        one_indexed = False

    # load the connectivity array: optional
    if 'face_face_connectivity' in names_mapping:
        face_face_connectivity = nc.variables[names_mapping['face_face_connectivity']]
        # does it need to be transposed?
        if face_face_connectivity.shape[0] <= face_face_connectivity.shape[1]:  # assume there are more than three triangles...
            # fortran order -- needs to be transposed
            face_face_connectivity = face_face_connectivity[:].T
        else:
            face_face_connectivity = face_face_connectivity[:]
        if one_indexed:
            face_face_connectivity -= 1
        ug.face_face_connectivity = face_face_connectivity[:]


    # load the center points of the faces: optional
    if 'face_coordinates_lon' in names_mapping and 'face_coordinates_lon' in names_mapping:
        ug.face_coordinates = np.zeros((len(ug.faces), 2), dtype=lon.dtype)
        ug.face_coordinates[:,0] = nc.variables[names_mapping['face_coordinates_lon']][:]
        ug.face_coordinates[:,1] = nc.variables[names_mapping['face_coordinates_lat']][:]


    print "checking bounds"
    if 'boundaries' in names_mapping: # optional
        print "loading boundaries"
        ## fixme --  this one is weird and non-conforming....
        boundaries = nc.variables[names_mapping['boundaries']][:,:2]
        if one_indexed:
            boundaries -= 1
        ug.boundaries = boundaries # ignoring the second two fields -- what are they???


    return ug


if __name__ == "__main__":

    names_mapping = {
                     'attribute_check': ('grid_type','triangular'),
                     'nodes_lon': 'lon',
                     'nodes_lat' : 'lat',
                     'faces' : 'nv',
                     'face_face_connectivity' : 'nbe',
                     'face_coordinates_lon' : 'lonc',
                     'face_coordinates_lat' : 'latc',
                     'boundaries' : 'bnd',
                    }


    ug = load_from_varnames("small_trigrid_example.nc", names_mapping)

    print ug





