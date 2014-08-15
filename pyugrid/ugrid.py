#!/usr/bin/env python

"""
ugrid classes

set of classes for working with unstructured model grids

The "ugrid" class is the base class: it stores everything in memory

It can read from and write to netcdf files in the UGRID format.

It may be able to reference a netcdf file at some point, rather than storing directly in memory.

NOTE: only support for triangular mesh grids at the moment

"""

import numpy as np

from . import read_netcdf
# used for simple locate_face test
#from py_geometry.cy_point_in_polygon import point_in_poly as point_in_tri
from .util import point_in_tri
from .data_set import DataSet

IND_DT = np.int32 ## datatype used for indexes -- might want to change for 64 bit some day.
NODE_DT = np.float64 ## datatype used for node coordinates




class UGrid(object):
    """
    a basic class to hold an unstructured grid (triangular mesh)

    the internal structure mirrors the netcdf data standard.
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
                 data = None,
                 mesh_name = "mesh",
                 ):
        """
        ugrid class -- holds, saves, etc. an unstructured grid

        :param nodes=None : the coordinates of the nodes -- (NX2) float array
        :param faces=None : the faces of the grid -- (NX3) integer array of indexes into the nodes array
        :param edges=None : the edges of the grid -- (NX2) integer array of indexes into the nodes array
        :param boundaries=None: specification of the boundaries -- are usually a subset of edges
                                where boundary condition information, etc is stored.
                                (NX2) integer array of indexes into the nodes array
        :type boundaries: numpy array of integer dtype


        :param face_face_connectivity=None: connectivity arrays
        :param face_edge_connectivity=None: connectivity arrays

        :param edge_coordinates=None: representative coordinate of the edges
        :param face_coordinates=None: representative coordinate of the faces
        :param boundary_coordinates=None: representative coordinate of the boundaries
        

        :param data = None: associated dataset
        :type data: dict of DataSet objects

        :param mesh_name = "mesh": optional name for the mesh 
        :type string: 

        often this is too much data to pass in as literals -- so usually
        specialized constructors will be used instead (load from file, etc.)
        """

        self.nodes = nodes
        self.faces = faces
        self.edges = edges
        self.boundaries = boundaries

        self.face_face_connectivity = face_face_connectivity
        self.face_edge_connectivity = face_edge_connectivity

        self.edge_coordinates=edge_coordinates
        self.face_coordinates=face_coordinates
        self.boundary_coordinates=boundary_coordinates


        self.mesh_name = mesh_name

        # the data associated with the grid
        # should be a dict of DataSet objects

        self._data = {} # the data associated with the grid
        if data is not None:
            for dataset in data.values():
                self.add_data(dataset)

    @classmethod
    def from_ncfile(klass, nc_url, mesh_name=None, load_data=False):
        """
        create a UGrid object from a netcdf file name (or opendap url)

        :param nc_url: the filename or OpenDap url you want to load

        :param mesh_name=None: the name of the mesh you want. If None, then
                               you'll get the only mesh in the file. If there
                               is more than one mesh in the file, a ValueError
                               Will be raised
        :param load_data=False: flag to indicate whether you want to load the associated
                                data or not. The mesh will be loaded in any case. If False,
                                only the mesh will be loaded. If True, then all the data
                                associated with the mesh will be loaded. This could be huge!
        :type load_data: boolean

        """
        grid = klass()
        read_netcdf.load_grid_from_ncfilename(nc_url, grid, mesh_name, load_data)
        return grid

    @classmethod
    def from_nc_dataset(klass, nc, mesh_name=None, load_data=False):
        """
        create a UGrid object from a netcdf file (or opendap url)

        :param nc: An already open Dataset object
        :type nc: netCDF4.DataSet

        :param mesh_name=None: the name of the mesh you want. If None, then
                               you'll get the only mesh in the file. If there
                               is more than one mesh in the file, a ValueError
                               Will be raised

        :param load_data=False: flag to indicate whether you want to load the associated
                                data or not. The mesh will be loaded in any case. If False,
                                only the mesh will be loaded. If True, then all the data
                                associated with the mesh will be loaded. This could be huge!

        :type load_data: boolean

        """
        grid = klass()
        read_netcdf.load_grid_from_nc_dataset(nc, grid, mesh_name, load_data)
        return grid

    def check_consistent(self):
        """
        check if the various data is consistent: the edges and faces reference
        existing nodes, etc.
        """
        raise NotImplimentedError
    
    @property
    def num_vertices(self):
        """
        number of vertices in a face
        """
        if self._faces is None:
            return None
        else:
            return self._faces.shape[1]

    @property
    def nodes(self):
        return self._nodes
    @nodes.setter
    def nodes(self, nodes_coords):
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
        if nodes_coords is None:
            self.nodes = np.zeros((0,2), dtype=NODE_DT)
        else:    
            self._nodes = np.asarray(nodes_coords, dtype=NODE_DT)

    @nodes.deleter
    def nodes(self):
        ## if there are no nodes, there can't be anything else
        self._nodes = np.zeros((0,2), dtype=NODE_DT)
        self._edges = None
        self._faces = None
        self._boundaries = None

    @property
    def faces(self):
        return self._faces
    @faces.setter
    def faces(self, faces_indexes):
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
        if faces_indexes is not None:
            self._faces = np.asarray(faces_indexes, dtype=IND_DT)
        else:
            self._faces = None
            # other things are no longer valid
            self._face_face_connectivity = None
            self._face_edge_connectivity = None
    @faces.deleter
    def faces(self):
        self._faces = None
        self._faces = None
        # other things are no longer valid
        self._face_face_connectivity = None
        self._face_edge_connectivity = None
        self.edge_coordinates = None

    @property
    def edges(self):
        return self._edges
    @edges.setter
    def edges(self, edges_indexes):
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
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
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
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
        ## add more checking?
        if face_face_connectivity is not None:
            face_face_connectivity = np.asarray(face_face_connectivity, dtype=IND_DT)
            if face_face_connectivity.shape != (len(self.faces), self.num_vertices):
                raise ValueError("face_face_connectivity must be size (num_faces, %i)"%self.num_vertices)
        self._face_face_connectivity = face_face_connectivity
    @face_face_connectivity.deleter
    def face_face_connectivity(self):
        self._face_face_connectivity = None

    @property
    def face_edge_connectivity(self):
        return self._face_edge_connectivity
    @face_face_connectivity.setter
    def face_edge_connectivity(self, face_edge_connectivity):
        ## add more checking?
        if face_edge_connectivity is not None:
            face_edge_connectivity = np.asarray(face_edge_connectivity, dtype=IND_DT)
            if face_edge_connectivity.shape != (len(self.faces), self.num_vertices):
                raise ValueError("face_face_connectivity must be size (num_face, %i)"%self.num_vertices)
        self._face_edge_connectivity = face_edge_connectivity
    @face_edge_connectivity.deleter
    def face_edge_connectivity(self):
        self._face_edge_connectivity = None


    @property
    def data(self):
        """
        dict of data associated with the data arrays
        You can't set this -- must use UGrid.add_data()
        """
        return self._data

    def add_data(self, data_set):
        """
        Add a dataset to the data dict

        :param data_set: the DataSet object to add. its name will be the key in the data dict.
        :type data_set: a ugrid.DataSet object

        some sanity checking is done to make sure array sizes are correct.

        """
        # do a size check:
        if data_set.location == 'node':
            if len(data_set.data) != len(self.nodes):
                raise ValueError("length of data array must match the number of nodes")
        elif data_set.location == 'edge':
            if len(data_set.data) != len(self.edges):
                raise ValueError("length of data array must match the number of edges")
        elif data_set.location == 'face':
            if len(data_set.data) != len(self.faces):
                raise ValueError("length of data array must match the number of faces")
        elif data_set.location == 'boundary':
            if len(data_set.data) != len(self.boundaries):
                raise ValueError("length of data array must match the number of boundaries")
        else:
            raise ValueError("I don't know how to add data associated with '%s'"%data_set.location)
        self._data[data_set.name] = data_set

    def find_data_sets(self, standard_name, location=None):
        """
        Find all :py:class:`DataSet`s that match the specified standard name
        
        :param str standard_name: the standard name attribute (based on the UGRID conventions)
        
        :keyword location: optional attribute location to narrow the returned
                           :py:class:`DataSet`s (one of 'node', 'edge', 'face', or 'boundary').
        
        :return: set of matching :py:class:`DataSet`s
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


    def locate_face_simple(self, point):
        """
        returns the index of the face that the point is in

        returns None if the point is not in the mesh

        : param point :  the point that you want to locate -- (x, y)

        this is a very simple, look through all the faces search.
        It is slow ( O(N) ), but should be robust
        """
        for i, face in enumerate(self._faces):
            f = self._nodes[face]
            if point_in_poly(f, point):
                return i
        return None

    def build_face_face_connectivity(self):
        """
        builds the face_face_connectivity array:
        essentially giving the neighbors of each triangle

        note: arbitrary order and CW vs CCW may not be consistent
        """        
        num_vertices = self.num_vertices
        num_faces = self.faces.shape[0]
        face_face = np.zeros( (num_faces, num_vertices), dtype=IND_DT )
        face_face += -1 # fill with -1

        # loop through all the triangles to find the matching edges:
        edges = {} # dict to store the edges in 
        for i, face in enumerate(self.faces):
            # loop through edges of the triangle:
            for j in range(num_vertices):
                if j < self.num_vertices-1:
                    edge = (face[j], face[j+1])
                else:
                    edge = (face[-1], face[0])
                if edge[0] > edge[1]: # sort the node numbers
                    edge = (edge[1], edge[0]) 
                # see if it is already in there
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i,j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j) # face num, edge_num
        self._face_face_connectivity = face_face

    def build_edges(self):
        """
        builds the edges array: all the edges defined by the triangles

        This will replace the existing edge array, if there is one.

        NOTE: arbitrary order -- should the order be preserved?
        """        
        num_vertices = self.num_vertices
        num_faces = self.faces.shape[0]
        face_face = np.zeros( (num_faces, num_vertices), dtype=IND_DT )
        face_face += -1 # fill with -1

        # loop through all the faces to find all the edges:
        edges = set() # use a set so no duplicates
        for i, face in enumerate(self.faces):
            # loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]: # flip them
                    edge = (edge[1], edge[0]) 
                edges.add(edge)
        self._edges = np.array(list(edges), dtype=IND_DT)

    def build_boundaries(self):
        """
        builds the boundary segments from the cell array

        It is assumed that -1 means no neighbor, which indicates a boundary
        
        This will over-write the existing boundaries array if there is one.

        This is a not-very-smart just loop through all the faces method.
        """
        boundaries = []
        for i, face in enumerate(self.face_face_connectivity):
            for j, neighbor in enumerate(face):
                if neighbor == -1:
                    if j == self.num_vertices-1:
                        bound = ( self.faces[i,-1], self.faces[i,0] )
                    else:
                        bound = ( self.faces[i,j], self.faces[i,j+1] )
                    boundaries.append(bound) 
        self.boundaries = boundaries

    def build_face_edge_connectivity(self):
        """
        Builds the face-edge connectivity array

        Not implemented yet.
        """
        raise NotImplimentedError

    def build_face_coordinates(self):
        """
        Builds the face_coordinates array, using the average of the
        nodes defining each face.

        Note that you may want a different definition of the face
        coordinates than this computes, but this is here to have
        an easy default.

        This will write-over an existing face_coordinates array

        Useful if you want this in the output file
        """
        face_coordinates = np.zeros( (len(self.faces),2), dtype=NODE_DT )
        ## fixme: there has got to be a way to vectorize this
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
        edge_coordinates = np.zeros( (len(self.edges),2), dtype=NODE_DT )
        ## fixme: there has got to be a way to vectorize this
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
        boundary_coordinates = np.zeros( (len(self.boundaries),2), dtype=NODE_DT )
        ## fixme: there has got to be a way to vectorize this
        for i, bound in enumerate(self.boundaries):
            coords = self.nodes[bound]
            boundary_coordinates[i] = coords.mean(axis=0)
        self.boundary_coordinates = boundary_coordinates


    def save_as_netcdf(self, filepath):
        """
        save the ugrid object as a netcdf file

        :param filepath: path to file you want o save to.
                         An existing one will be clobbered if it already exists.

        follows the convernsion established by the netcdf UGRID working group:

        http://publicwiki.deltares.nl/display/NETCDF/Deltares+CF+proposal+for+Unstructured+Grid+data+model

        """
        mesh_name = self.mesh_name


        from netCDF4 import Dataset as ncDataset
        from netCDF4 import num2date, date2num
        # create a new netcdf file
        with ncDataset(filepath, mode="w", clobber=True) as nclocal:

            nclocal.createDimension(mesh_name+'_num_node', len(self.nodes) )
            if self._edges is not None:
                nclocal.createDimension(mesh_name+'_num_edge', len(self.edges) )
            if self._boundaries is not None:
                nclocal.createDimension(mesh_name+'_num_boundary', len(self.boundaries) )
            if self._faces is not None:
                nclocal.createDimension(mesh_name+'_num_face', len(self.faces) )
                nclocal.createDimension(mesh_name+'_num_vertices', self.faces.shape[1] )
            nclocal.createDimension('two', 2)

            #mesh topology
            mesh = nclocal.createVariable(mesh_name, IND_DT, (), )
            mesh.cf_role = "mesh_topology" 
            mesh.long_name = "Topology data of 2D unstructured mesh" 
            mesh.topology_dimension = 2 
            mesh.node_coordinates = "{0}_node_lon {0}_node_lat".format(mesh_name) 

            if self.edges is not None:
                mesh.edge_node_connectivity = mesh_name+"_edge_nodes"  ## attribute required if variables will be defined on edges
                if self.edge_coordinates is not None:
                    mesh.edge_coordinates =   "{0}_edge_lon {0}_edge_lat".format(mesh_name)  ## optional attribute (requires edge_node_connectivity)
            if self.faces is not None:
                mesh.face_node_connectivity = mesh_name+"_face_nodes"
                if self.face_coordinates is not None: 
                    mesh.face_coordinates = "{0}_face_lon {0}_face_lat".format(mesh_name) ##  optional attribute
            if self.face_edge_connectivity is not None:
                mesh.face_edge_connectivity = mesh_name+"_face_edges"  ## optional attribute (requires edge_node_connectivity)
            if self.face_face_connectivity is not None:
                mesh.face_face_connectivity = mesh_name+"_face_links"  ## optional attribute
            if self.boundaries is not None:
                mesh.boundary_node_connectivity = mesh_name+"_boundary_nodes"

            ## fixme: This could be re-factored to be more generic, rather than separate for each type of data
            ##        see the coordinates example below
            if self.faces is not None:
                face_nodes = nclocal.createVariable(mesh_name+"_face_nodes",
                                                    IND_DT,
                                                    (mesh_name+'_num_face', mesh_name+'_num_vertices'),
                                                    )
                face_nodes[:] = self.faces

                face_nodes.cf_role = "face_node_connectivity"
                face_nodes.long_name = "Maps every triangular face to its three corner nodes."
                face_nodes.start_index = 0 ;

            if self.edges is not None:
                edge_nodes = nclocal.createVariable(mesh_name+"_edge_nodes",
                                                    IND_DT,
                                                    (mesh_name+'_num_edge', 'two'),
                                                    )
                edge_nodes[:] = self.edges

                edge_nodes.cf_role = "edge_node_connectivity"
                edge_nodes.long_name = "Maps every edge to the two nodes that it connects."
                edge_nodes.start_index = 0 ;

            if self.boundaries is not None:
                boundary_nodes = nclocal.createVariable(mesh_name+"_boundary_nodes",
                                                        IND_DT,
                                                        (mesh_name+'_num_boundary', 'two'),
                                                        )
                boundary_nodes[:] = self.boundaries

                boundary_nodes.cf_role = "boundary_node_connectivity"
                boundary_nodes.long_name = "Maps every boundary segment to the two nodes that it connects."
                boundary_nodes.start_index = 0 ;

            ## optional "coordinate variables"
            for location in ['face', 'edge', 'boundary']:
                if getattr(self, "{0}_coordinates".format(location)) is not None:
                    for axis, ind in [('lat',1), ('lon',0)]:
                        var = nclocal.createVariable("{0}_{1}_{2}".format(mesh_name, location, axis),
                                                     NODE_DT,
                                                     dimensions=("{0}_num_{1}".format(mesh_name, location)),
                                                    )
                        var[:] = getattr(self, "{0}_coordinates".format(location))[:,ind]
                        ## attributes of the variable
                        var.standard_name = "longitude" if axis == 'lon' else 'latitude'
                        var.units = "degrees_east" if axis == 'lon' else 'degrees_north'
                        var.long_name = "Characteristics {0} of 2D mesh {1}".format(var.standard_name, location)

            ## the node data
            node_lon = nclocal.createVariable(mesh_name+'_node_lon',
                                              self._nodes.dtype,
                                              (mesh_name+'_num_node',),
                                              chunksizes=(len(self.nodes), ),
                                              #zlib=False,
                                              #complevel=0,
                                              )
            node_lon[:] = self.nodes[:,0]
            node_lon.standard_name = "longitude"
            node_lon.long_name = "Longitude of 2D mesh nodes."
            node_lon.units = "degrees_east"

            node_lat = nclocal.createVariable(mesh_name+'_node_lat',
                                              self._nodes.dtype,
                                              (mesh_name+'_num_node',),
                                              chunksizes=(len(self.nodes), ),
                                              #zlib=False,
                                              #complevel=0,
                                              )
            node_lat[:] = self.nodes[:,1]
            node_lat.standard_name = "latitude"
            node_lat.long_name = "Latitude of 2D mesh nodes."
            node_lat.units = "degrees_north"

            ## write the associated data
            for dataset in self.data.values():
                if dataset.location == 'node':
                    shape = (mesh_name + '_num_node',)
                    coordinates = "{0}_node_lon {0}_node_lat".format(mesh_name)
                    chunksizes = (len(self.nodes), )
                elif dataset.location == 'face':
                    shape = (mesh_name + '_num_face',)
                    coordinates = "{0}_face_lon {0}_face_lat".format(mesh_name) if self.face_coordinates is not None else None
                    chunksizes = (len(self.faces), )
                elif dataset.location == 'edge':
                    shape = (mesh_name + '_num_edge',)
                    coordinates = "{0}_edge_lon {0}_edge_lat".format(mesh_name) if self.edge_coordinates is not None else None
                    chunksizes = (len(self.edges), )
                elif dataset.location == 'boundary':
                    shape = (mesh_name + '_num_boundary',)
                    coordinates = "{0}_boundary_lon {0}_boundary_lat".format(mesh_name) if self.boundary_coordinates is not None else None
                    chunksizes = (len(self.boundaries), )
                data_var = nclocal.createVariable(dataset.name,
                                                  dataset.data.dtype,
                                                  shape,
                                                  chunksizes=chunksizes,
                                                  #zlib=False,
                                                  #complevel=0,
                                                  )
                data_var[:] = dataset.data
                ## add the standard attributes:
                data_var.location = dataset.location
                data_var.mesh = mesh_name
                if coordinates is not None:
                    data_var.coordinates = coordinates
                ## add the extra attributes
                for att_name, att_value in dataset.attributes.items():
                    setattr(data_var, att_name, att_value)

            nclocal.sync()









