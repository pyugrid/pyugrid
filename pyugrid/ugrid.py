#!/usr/bin/env python

"""
UGrid base class.

NOTE: Supports only triangular mesh grids at the moment.

"""

from __future__ import (absolute_import, division, print_function)

import numpy as np
from netCDF4 import Dataset

from .read_netcdf import find_mesh_names, is_valid_mesh, grid_defs, coord_defs
from .util import point_in_tri
from .uvar import UVar


__all__ = ['UGrid',
           'load_grid']

# Data type used for indexes. Might want to change for 64 bit some day.
IND_DT = np.int32
NODE_DT = np.float64  # Data type used for node coordinates.


class UGrid(object):
    """
    A basic class to hold an unstructured grid (triangular mesh).

    The internal structure mirrors the netcdf data standard.
    """

    def __init__(self,
                 nodes=None,
                 faces=None,
                 edges=None,
                 boundaries=None,
                 face_face_connectivity=None,
                 face_edge_connectivity=None,
                 edge_coordinates=None,
                 face_coordinates=None,
                 boundary_coordinates=None,
                 data=None,
                 mesh_name="mesh",
                 ):
        """
        ugrid class -- holds, saves, etc. an unstructured grid

        :param nodes=None : the coordinates of the nodes
        :type nodes: (NX2) array of floats

        :param faces=None : the faces of the grid. Indexes for the nodes array.
        :type faces: (NX3) array of integers

        :param edges=None : the edges of the grid. Indexes for the nodes array.
        :type edges: (NX2) array of integers

        :param boundaries=None: specification of the boundaries are usually a
                                subset of edges where boundary condition
                                information, etc is stored.
                                (NX2) integer array of indexes for the nodes
                                array.
        :type boundaries: numpy array of integers

        :param face_face_connectivity=None: connectivity arrays.
        :param face_edge_connectivity=None: connectivity arrays.

        :param edge_coordinates=None: representative coordinate of the edges.
        :param face_coordinates=None: representative coordinate of the faces.
        :param boundary_coordinates=None: representative coordinate of the
                                          boundaries.

        :param edge_coordinates=None: representative coordinate of the edges
        :type edge_coordinates: (NX2) array of floats

        :param face_coordinates=None: representative coordinate of the faces
                                      (NX2) float array
        :type face_coordinates: (NX2) array of floats


        :param boundary_coordinates=None: representative coordinate of the
                                          boundaries
        :type boundary_coordinates: (NX2) array of floats


        :param data = None: associated variables
        :type data: dict of UVar objects

        :param mesh_name = "mesh": optional name for the mesh
        :type mesh_name: string

        Often this is too much data to pass in as literals -- so usually
        specialized constructors will be used instead (load from file, etc).

        """
        self.nodes = nodes
        self.faces = faces
        self.edges = edges
        self.boundaries = boundaries

        self.face_face_connectivity = face_face_connectivity
        self.face_edge_connectivity = face_edge_connectivity

        self.edge_coordinates = edge_coordinates
        self.face_coordinates = face_coordinates
        self.boundary_coordinates = boundary_coordinates

        self.mesh_name = mesh_name

        # The data associated with the grid should be a dict of UVar objects.
        self._data = {}
        if data is not None:
            for dataset in data.values():
                self.add_data(dataset)

        self._kdtree = None

    @classmethod
    def load_grid(klass, nc, mesh_name=None, load_data=False):
        """
        Loads UGrid object from a `netCDF4.Dataset` object or file name/URL,
        adding the data to the passed-in grid object.

        It will load the mesh specified, or look for the first one it finds if
        none is specified.

        :param nc: The grid data to be loaded up.
        :type nc: `netCDF4.Dataset` object or netCDF/URL string

        :param mesh_name=None: name of the mesh to load.
        :type mesh_name: string

        :param load_data=False: flag to indicate whether you want to load the
                                associated data or not. The mesh will be loaded
                                in any case. If False, only the mesh will be
                                loaded. If True, then all the data associated
                                with the mesh will be loaded.
                                This could be huge!
        :type load_data: boolean

        """
        ugrid = klass()

        if isinstance(nc, Dataset):
            pass
        else:
            nc = Dataset(nc, 'r')

        # FIXME: factor this code out in a compliance-checker function/module.
        # Get the mesh_name.
        if mesh_name is None:
            # Find the mesh.
            meshes = find_mesh_names(nc)
            if len(meshes) == 0:
                msg = "There are no standard-conforming meshes in {}".format
                # FIXME: The use of `nc.filepath`
                # fails in a non-cythonized install of netCDF4-python!
                raise ValueError(msg(nc.filepath))
            if len(meshes) > 1:
                msg = "There is more than one mesh in the file: {!r}".format
                raise ValueError(msg(meshes))
            mesh_name = meshes[0]
        else:
            if not is_valid_mesh(nc, mesh_name):
                msg = "Mesh: {} is not in {}".format
                raise ValueError(msg(mesh_name, nc.filepath))

        ugrid.mesh_name = mesh_name
        mesh_var = nc[mesh_name]

        # Load the coordinate variables.
        for defs in coord_defs:
            try:
                coord_names = mesh_var.getncattr(defs['role']).strip().split()
                coord_vars = [nc.variables[name] for name in coord_names]
            except AttributeError:
                if defs['required']:
                    msg = "Mesh variable must include {} attribute.".format
                    raise ValueError(msg(defs['role']))
                continue
            except KeyError:
                msg = ("File must include {} variables for {} "
                       "named in mesh variable.").format
                raise ValueError(msg(coord_names, defs['role']))

            coord_vars = [nc.variables[name] for name in coord_names]
            num_node = len(coord_vars[0])
            nodes = np.empty((num_node, 2), dtype=np.float64)
            for var in coord_vars:
                # CF does not require a standard_name.
                try:
                    standard_name = var.standard_name
                except AttributeError:
                    # Look for units.
                    try:
                        units = var.units
                    except AttributeError:
                        msg = ("The {} variable doesn't contain units "
                               "attribute: required by CF").format
                        raise ValueError(msg(var))
                    # CF accepted units attributes for longitude.
                    if units in ('degrees_east', 'degree_east', 'degree_E',
                                 'degrees_E', 'degreeE', 'degreesE'):
                            standard_name = 'longitude'
                    # CF accepted units attributes for longitude.
                    elif units in ('degrees_north', 'degree_north', 'degree_N',
                                   'degrees_N', 'degreeN', 'degreesN'):
                            standard_name = 'latitude'
                    else:
                        msg = ("{} variable's units value ({}) doesn't look "
                               "like latitude or longitude").format
                        raise ValueError(msg(var, units))
                # FIXME: `nodes_lon` and `nodes_lat`? We will be compatible
                # with other libraries and we can just pass the netCDF4-python
                # var without the slicing.
                # (Getting a deferred lazy load for free.)
                if standard_name == 'latitude':
                    nodes[:, 1] = var[:]
                elif standard_name == 'longitude':
                    nodes[:, 0] = var[:]
                else:
                    raise ValueError('Node coordinates standard_name is '
                                     'neither "longitude" nor "latitude".')
            setattr(ugrid, defs['grid_attr'], nodes)

        # Load assorted connectivity arrays.
        for defs in grid_defs:
            try:
                try:
                    var = nc.variables[mesh_var.getncattr(defs['role'])]
                except AttributeError:  # This connectivity array isn't there.
                    continue
                array = var[:, :]
                # Fortran order, instead of C order, transpose the array
                # logic below will fail for 3 node or two edge grids.
                if array.shape[0] == defs['num_ind']:
                    array = array.T
                try:
                    start_index = int(var.start_index)
                except AttributeError:
                    start_index = 0
                if start_index >= 1:
                    array -= start_index
                    # Check for flag value.
                    try:
                        # FIXME: This won't work for more than one flag value.
                        flag_value = var.flag_values
                        array[array == flag_value-start_index] = flag_value
                    except AttributeError:
                        pass
                setattr(ugrid, defs['grid_attr'], array)
            except KeyError:
                pass  # OK not to have this...

        # Load the associated data:
        if load_data:
            # Look for data arrays -- they should have a "location" attribute.
            for name, var in nc.variables.items():
                # Data Arrays should have "location" and "mesh" attributes.
                try:
                    location = var.location
                    # The mesh attribute should match the mesh we're loading:
                    if var.mesh != mesh_name:
                        continue
                except AttributeError:
                    continue

                # Get the attributes.
                # FIXME: Is there a way to get the attributes directly?
                attributes = {n: var.getncattr(n) for n in var.ncattrs()
                              if n not in ('location', 'coordinates', 'mesh')}

                # Trick with the name: FIXME: Is this a good idea?
                name = name.lstrip(mesh_name).lstrip('_')
                uvar = UVar(name, data=var[:],
                            location=location, attributes=attributes)
                ugrid.add_data(uvar)
        return ugrid

    def check_consistent(self):
        """
        Check if the various data is consistent: the edges and faces reference
        existing nodes, etc.

        """
        raise NotImplementedError

    @property
    def num_vertices(self):
        """Number of vertices in a face."""
        if self._faces is None:
            return None
        else:
            return self._faces.shape[1]

    @property
    def nodes(self):
        return self._nodes

    @nodes.setter
    def nodes(self, nodes_coords):
        if nodes_coords is None:
            self.nodes = np.zeros((0, 2), dtype=NODE_DT)
        else:
            self._nodes = np.asarray(nodes_coords, dtype=NODE_DT)

    @nodes.deleter
    def nodes(self):
        self._nodes = np.zeros((0, 2), dtype=NODE_DT)
        self._edges = None
        self._faces = None
        self._boundaries = None

    @property
    def faces(self):
        return self._faces

    @faces.setter
    def faces(self, faces_indexes):
        if faces_indexes is not None:
            self._faces = np.asarray(faces_indexes, dtype=IND_DT)
        else:
            self._faces = None
            # Other things are no longer valid.
            self._face_face_connectivity = None
            self._face_edge_connectivity = None

    @faces.deleter
    def faces(self):
        self._faces = None
        self._faces = None
        # Other things are no longer valid.
        self._face_face_connectivity = None
        self._face_edge_connectivity = None
        self.edge_coordinates = None

    @property
    def edges(self):
        return self._edges

    @edges.setter
    def edges(self, edges_indexes):
        if edges_indexes is not None:
            self._edges = np.asarray(edges_indexes, dtype=IND_DT)
        else:
            self._edges = None
            self._face_edge_connectivity = None

    @edges.deleter
    def edges(self):
        self._edges = None
        self._face_edge_connectivity = None
        self.edge_coordinates = None

    @property
    def boundaries(self):
        return self._boundaries

    @boundaries.setter
    def boundaries(self, boundaries_indexes):
        if boundaries_indexes is not None:
            self._boundaries = np.asarray(boundaries_indexes, dtype=IND_DT)
        else:
            self._boundaries = None

    @boundaries.deleter
    def boundaries(self):
        self._boundaries = None
        self.boundary_coordinates = None

    @property
    def face_face_connectivity(self):
        return self._face_face_connectivity

    @face_face_connectivity.setter
    def face_face_connectivity(self, face_face_connectivity):
        """TODO: Write docstrings."""
        # Add more checking?
        if face_face_connectivity is not None:
            face_face_connectivity = np.asarray(face_face_connectivity,
                                                dtype=IND_DT)
            if face_face_connectivity.shape != (len(self.faces),
                                                self.num_vertices):
                msg = ("face_face_connectivity must be size "
                       "(num_faces, {})").format
                raise ValueError(msg(self.num_vertices))
        self._face_face_connectivity = face_face_connectivity

    @face_face_connectivity.deleter
    def face_face_connectivity(self):
        self._face_face_connectivity = None

    @property
    def face_edge_connectivity(self):
        return self._face_edge_connectivity

    @face_edge_connectivity.setter
    def face_edge_connectivity(self, face_edge_connectivity):
        """TODO: Write docstrings."""
        # Add more checking?
        if face_edge_connectivity is not None:
            face_edge_connectivity = np.asarray(face_edge_connectivity,
                                                dtype=IND_DT)
            if face_edge_connectivity.shape != (len(self.faces),
                                                self.num_vertices):
                msg = ("face_face_connectivity must be size "
                       "(num_face, {})").format
                raise ValueError(msg(self.num_vertices))
        self._face_edge_connectivity = face_edge_connectivity

    @face_edge_connectivity.deleter
    def face_edge_connectivity(self):
        self._face_edge_connectivity = None

    @property
    def data(self):
        """
        dict of data associated with the data arrays
        You can't set this -- must use UGrid.add_data().

        """
        return self._data

    def add_data(self, uvar):
        """
        Add a UVar to the data dict

        :param uvar: the UVar object to add.
                     Its name will be the key in the data dict.
        :type uvar: a ugrid.UVar object

        Some sanity checking is done to make sure array sizes are correct.

        """
        # Size check:
        if uvar.location == 'node':
            if len(uvar.data) != len(self.nodes):
                raise ValueError("length of data array must match "
                                 "the number of nodes")
        elif uvar.location == 'edge':
            if len(uvar.data) != len(self.edges):
                raise ValueError("length of data array must match "
                                 "the number of edges")
        elif uvar.location == 'face':
            if len(uvar.data) != len(self.faces):
                raise ValueError("length of data array must match "
                                 "the number of faces")
        elif uvar.location == 'boundary':
            if len(uvar.data) != len(self.boundaries):
                raise ValueError("length of data array must match "
                                 "the number of boundaries")
        else:
            msg = "Can't add data associated with '{}'".format
            raise ValueError(msg(uvar.location))
        self._data[uvar.name] = uvar

    def find_uvars(self, standard_name, location=None):
        """
        Find all :py:class:`UVar`s that match the specified standard name

        :param str standard_name: the standard name attribute.
                                  Based on the UGRID conventions.

        :keyword location: optional attribute location to narrow the returned
        :py:class : `UVar`s (one of 'node', 'edge', 'face', or 'boundary').

        :return: set of matching :py:class:`UVar`s

        """
        found = set()
        for ds in self._data.values():
            if not ds.attributes or 'standard_name' not in ds.attributes:
                continue
            if ds.attributes['standard_name'] == standard_name:
                if location is not None and ds.location != location:
                    continue
                found.add(ds)
        return found

    def locate_nodes(self, points):
        """
        Returns the index of the closest nodes to the input locations.

        :param points: the lons/lats of locations you want the nodes
                       closest to.
        :type point: a (N, 2) ndarray of points
                     (or something that can be converted).

        :returns: the index of the closest node.

        """
        if self._kdtree is None:
            self._build_kdtree()

        node_inds = self._kdtree.query(points, k=1)[1]
        return node_inds

    def _build_kdtree(self):
        from scipy.spatial import cKDTree
        self._kdtree = cKDTree(self.nodes)

    def locate_face_simple(self, point):
        """
        Returns the index of the face that the point is in.

        Returns `None` if the point is not in the mesh.

        :param point:  the point that you want to locate -- (x, y).

        This is a very simple, look through all the faces search.
        It is slow ( O(N) ), but should be robust

        """
        for i, face in enumerate(self._faces):
            f = self._nodes[face]
            if point_in_tri(f, point):
                return i
        return None

    def build_face_face_connectivity(self):
        """
        Builds the face_face_connectivity array:
        essentially giving the neighbors of each triangle.

        Note: arbitrary order and CW vs CCW may not be consistent.

        """
        num_vertices = self.num_vertices
        num_faces = self.faces.shape[0]
        face_face = np.zeros((num_faces, num_vertices), dtype=IND_DT)
        face_face += -1  # Fill with -1.

        # Loop through all the triangles to find the matching edges:
        edges = {}  # dict to store the edges.
        for i, face in enumerate(self.faces):
            # Loop through edges of the triangle:
            for j in range(num_vertices):
                if j < self.num_vertices-1:
                    edge = (face[j], face[j+1])
                else:
                    edge = (face[-1], face[0])
                if edge[0] > edge[1]:  # Sort the node numbers.
                    edge = (edge[1], edge[0])
                # See if it is already in there.
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i, j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j)  # face num, edge_num.
        self._face_face_connectivity = face_face

    def build_edges(self):
        """
        Builds the edges array: all the edges defined by the triangles

        This will replace the existing edge array, if there is one.

        NOTE: arbitrary order -- should the order be preserved?

        """
        num_vertices = self.num_vertices
        num_faces = self.faces.shape[0]
        face_face = np.zeros((num_faces, num_vertices), dtype=IND_DT)
        face_face += -1  # Fill with -1.

        # Loop through all the faces to find all the edges:
        edges = set()  # Use a set so no duplicates.
        for i, face in enumerate(self.faces):
            # Loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]:  # Flip them.
                    edge = (edge[1], edge[0])
                edges.add(edge)
        self._edges = np.array(list(edges), dtype=IND_DT)

    def build_boundaries(self):
        """
        Builds the boundary segments from the cell array.

        It is assumed that -1 means no neighbor, which indicates a boundary

        This will over-write the existing boundaries array if there is one.

        This is a not-very-smart just loop through all the faces method.

        """
        boundaries = []
        for i, face in enumerate(self.face_face_connectivity):
            for j, neighbor in enumerate(face):
                if neighbor == -1:
                    if j == self.num_vertices-1:
                        bound = (self.faces[i, -1], self.faces[i, 0])
                    else:
                        bound = (self.faces[i, j], self.faces[i, j+1])
                    boundaries.append(bound)
        self.boundaries = boundaries

    def build_face_edge_connectivity(self):
        """
        Builds the face-edge connectivity array

        Not implemented yet.

        """
        raise NotImplementedError

    def build_face_coordinates(self):
        """
        Builds the face_coordinates array, using the average of the
        nodes defining each face.

        Note that you may want a different definition of the face
        coordinates than this computes, but this is here to have
        an easy default.

        This will write-over an existing face_coordinates array.

        Useful if you want this in the output file.

        """
        face_coordinates = np.zeros((len(self.faces), 2), dtype=NODE_DT)
        # FIXME: there has got to be a way to vectorize this.
        for i, face in enumerate(self.faces):
            coords = self.nodes[face]
            face_coordinates[i] = coords.mean(axis=0)
        self.face_coordinates = face_coordinates

    def build_edge_coordinates(self):
        """
        Builds the face_coordinates array, using the average of the
        nodes defining each edge.

        Note that you may want a different definition of the edge
        coordinates than this computes, but this is here to have
        an easy default.


        This will write-over an existing face_coordinates array

        Useful if you want this in the output file

        """
        edge_coordinates = np.zeros((len(self.edges), 2), dtype=NODE_DT)
        # FIXME: there has got to be a way to vectorize this.
        for i, edge in enumerate(self.edges):
            coords = self.nodes[edge]
            edge_coordinates[i] = coords.mean(axis=0)
        self.edge_coordinates = edge_coordinates

    def build_boundary_coordinates(self):
        """
        Builds the boundary_coordinates array, using the average of the
        nodes defining each boundary segment.

        Note that you may want a different definition of the boundary
        coordinates than this computes, but this is here to have
        an easy default.

        This will write-over an existing face_coordinates array

        Useful if you want this in the output file

        """
        boundary_coordinates = np.zeros((len(self.boundaries), 2),
                                        dtype=NODE_DT)
        # FXIME: there has got to be a way to vectorize this.
        for i, bound in enumerate(self.boundaries):
            coords = self.nodes[bound]
            boundary_coordinates[i] = coords.mean(axis=0)
        self.boundary_coordinates = boundary_coordinates

    def save_as_netcdf(self, filepath):
        """
        Save the ugrid object as a netcdf file.

        :param filepath: path to file you want o save to.  An existing one
                         will be clobbered if it already exists.

        Follows the convention established by the netcdf UGRID working group:

        http://publicwiki.deltares.nl/display/NETCDF/Deltares+CF+proposal+for+Unstructured+Grid+data+model

        """
        mesh_name = self.mesh_name

        # Create a new netcdf file.
        with Dataset(filepath, mode='w', clobber=True) as nclocal:
            nclocal.createDimension(mesh_name+'_num_node', len(self.nodes))
            if self._edges is not None:
                nclocal.createDimension(mesh_name+'_num_edge', len(self.edges))
            if self._boundaries is not None:
                nclocal.createDimension(mesh_name+'_num_boundary',
                                        len(self.boundaries))
            if self._faces is not None:
                nclocal.createDimension(mesh_name+'_num_face', len(self.faces))
                nclocal.createDimension(mesh_name+'_num_vertices',
                                        self.faces.shape[1])
            nclocal.createDimension('two', 2)

            # Mesh topology.
            mesh = nclocal.createVariable(mesh_name, IND_DT, (),)
            mesh.cf_role = "mesh_topology"
            mesh.long_name = "Topology data of 2D unstructured mesh"
            mesh.topology_dimension = 2
            mesh.node_coordinates = "{0}_node_lon {0}_node_lat".format(mesh_name)  # noqa

            if self.edges is not None:
                # Attribute required if variables will be defined on edges.
                mesh.edge_node_connectivity = mesh_name + "_edge_nodes"
                if self.edge_coordinates is not None:
                    # Optional attribute (requires edge_node_connectivity).
                    coord = "{0}_edge_lon {0}_edge_lat".format
                    mesh.edge_coordinates = coord(mesh_name)
            if self.faces is not None:
                mesh.face_node_connectivity = mesh_name+"_face_nodes"
                if self.face_coordinates is not None:
                    # Optional attribute.
                    coord = "{0}_face_lon {0}_face_lat".format
                    mesh.face_coordinates = coord(mesh_name)
            if self.face_edge_connectivity is not None:
                # Optional attribute (requires edge_node_connectivity).
                mesh.face_edge_connectivity = mesh_name + "_face_edges"
            if self.face_face_connectivity is not None:
                # Optional attribute.
                mesh.face_face_connectivity = mesh_name + "_face_links"
            if self.boundaries is not None:
                mesh.boundary_node_connectivity = mesh_name+"_boundary_nodes"

            # FIXME: This could be re-factored to be more generic, rather than
            # separate for each type of data see the coordinates example below.
            if self.faces is not None:
                nc_create_var = nclocal.createVariable
                face_nodes = nc_create_var(mesh_name + "_face_nodes", IND_DT,
                                           (mesh_name + '_num_face',
                                            mesh_name + '_num_vertices'),)
                face_nodes[:] = self.faces

                face_nodes.cf_role = "face_node_connectivity"
                face_nodes.long_name = ("Maps every triangular face to "
                                        "its three corner nodes.")
                face_nodes.start_index = 0

            if self.edges is not None:
                nc_create_var = nclocal.createVariable
                edge_nodes = nc_create_var(mesh_name + "_edge_nodes", IND_DT,
                                           (mesh_name + '_num_edge', 'two'),)
                edge_nodes[:] = self.edges

                edge_nodes.cf_role = "edge_node_connectivity"
                edge_nodes.long_name = ("Maps every edge to the two "
                                        "nodes that it connects.")
                edge_nodes.start_index = 0

            if self.boundaries is not None:
                nc_create_var = nclocal.createVariable
                boundary_nodes = nc_create_var(mesh_name + "_boundary_nodes",
                                               IND_DT,
                                               (mesh_name + '_num_boundary',
                                                'two'),)
                boundary_nodes[:] = self.boundaries

                boundary_nodes.cf_role = "boundary_node_connectivity"
                boundary_nodes.long_name = ("Maps every boundary segment to "
                                            "the two nodes that it connects.")
                boundary_nodes.start_index = 0

            # Optional "coordinate variables."
            for location in ['face', 'edge', 'boundary']:
                loc = "{0}_coordinates".format(location)
                if getattr(self, loc) is not None:
                    for axis, ind in [('lat', 1), ('lon', 0)]:
                        nc_create_var = nclocal.createVariable
                        name = "{0}_{1}_{2}".format(mesh_name, location, axis)
                        dimensions = "{0}_num_{1}".format(mesh_name, location)
                        var = nc_create_var(name, NODE_DT,
                                            dimensions=(dimensions),)
                        loc = "{0}_coordinates".format(location)
                        var[:] = getattr(self, loc)[:, ind]
                        # Attributes of the variable.
                        var.standard_name = ("longitude" if axis == 'lon'
                                             else 'latitude')
                        var.units = ("degrees_east" if axis == 'lon'
                                     else 'degrees_north')
                        name = "Characteristics {0} of 2D mesh {1}".format
                        var.long_name = name(var.standard_name, location)

            # The node data.
            node_lon = nclocal.createVariable(mesh_name+'_node_lon',
                                              self._nodes.dtype,
                                              (mesh_name+'_num_node',),
                                              chunksizes=(len(self.nodes), ),
                                              # zlib=False,
                                              # complevel=0,
                                              )
            node_lon[:] = self.nodes[:, 0]
            node_lon.standard_name = "longitude"
            node_lon.long_name = "Longitude of 2D mesh nodes."
            node_lon.units = "degrees_east"

            node_lat = nclocal.createVariable(mesh_name+'_node_lat',
                                              self._nodes.dtype,
                                              (mesh_name+'_num_node',),
                                              chunksizes=(len(self.nodes), ),
                                              # zlib=False,
                                              # complevel=0,
                                              )
            node_lat[:] = self.nodes[:, 1]
            node_lat.standard_name = "latitude"
            node_lat.long_name = "Latitude of 2D mesh nodes."
            node_lat.units = "degrees_north"

            # Write the associated data.
            for dataset in self.data.values():
                if dataset.location == 'node':
                    shape = (mesh_name + '_num_node',)
                    coordinates = "{0}_node_lon {0}_node_lat".format(mesh_name)
                    chunksizes = (len(self.nodes), )
                elif dataset.location == 'face':
                    shape = (mesh_name + '_num_face',)
                    coord = "{0}_face_lon {0}_face_lat".format
                    coordinates = (coord(mesh_name) if self.face_coordinates
                                   is not None else None)
                    chunksizes = (len(self.faces), )
                elif dataset.location == 'edge':
                    shape = (mesh_name + '_num_edge',)
                    coord = "{0}_edge_lon {0}_edge_lat".format
                    coordinates = (coord(mesh_name) if self.edge_coordinates
                                   is not None else None)
                    chunksizes = (len(self.edges), )
                elif dataset.location == 'boundary':
                    shape = (mesh_name + '_num_boundary',)
                    coord = "{0}_boundary_lon {0}_boundary_lat".format
                    bcoord = self.boundary_coordinates
                    coordinates = (coord(mesh_name) if bcoord
                                   is not None else None)
                    chunksizes = (len(self.boundaries), )
                data_var = nclocal.createVariable(dataset.name,
                                                  dataset.data.dtype,
                                                  shape,
                                                  chunksizes=chunksizes,
                                                  # zlib=False,
                                                  # complevel=0,
                                                  )
                data_var[:] = dataset.data
                # Add the standard attributes:
                data_var.location = dataset.location
                data_var.mesh = mesh_name
                if coordinates is not None:
                    data_var.coordinates = coordinates
                # Add the extra attributes.
                for att_name, att_value in dataset.attributes.items():
                    setattr(data_var, att_name, att_value)
            nclocal.sync()


def load_grid(nc, mesh_name=None, load_data=False):
    """
    Loads UGrid object from a `netCDF4.Dataset` object or file name/URL,
    adding the data to the passed-in grid object.

    It will load the mesh specified, or look for the first one it finds if
    none is specified.

    :param nc: The grid data to be loaded up.
    :type nc: `netCDF4.Dataset` object or netCDF/URL string

    :param mesh_name=None: name of the mesh to load.
    :type mesh_name: string

    :param load_data=False: flag to indicate whether you want to load the
                            associated data or not. The mesh will be loaded
                            in any case. If False, only the mesh will be
                            loaded. If True, then all the data associated
                            with the mesh will be loaded.
                            This could be huge!
    :type load_data: boolean

    """
    if isinstance(nc, Dataset):
        pass
    else:
        nc = Dataset(nc, 'r')

    return UGrid.load_grid(nc, mesh_name=None, load_data=True)
