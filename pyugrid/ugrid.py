#!/usr/bin/env python

"""
ugrid classes

set of classes for working with unstructured model grids

The "ugrid" class is the base class: it stores eveything in memory

subclasses include a nc_ugrid, which points to a netcdf file (Or opendap url).

It provides the same API, but does not store the data in memory, rather
reading it on demand

NOTE: only full support for triangular mesh grids at the moment

"""

import numpy as np

# used for simple locate_face test
#from py_geometry.cy_point_in_polygon import point_in_poly as point_in_tri
from .util import point_in_tri

IND_DT = np.int32 ## datatype used for indexes -- might want to change for 64 bit some day.
NODE_DT = np.float64 ## datatype used for node coordinates


class data_set(object):
    """
    a class to hold the data assocated with nodes, edges, etc

    """
    ##fixme: is this neccessary? -- only if there is more in it later..
    def __init__(self, name, type='node', data=None):
        """
        create a data_set object
        name: the name of the data (depth, u_velocity, etc)
        type: the type of grid element: node, edge, or face the data is assigned to
        data: the data
        """
        self.name = name
        if type not in ['node', 'edge', 'face']:
            raise ValueError("type must be one of: 'node', 'edge', 'face'")
        self.type = type # must be 'node', 'edge', of 'face' (eventually 'volume')

        if data is None:
            self.data = np.zeros((0,), dtype=np.float) # could be any data type
        else:
            self.data = data

class data_set_indexed(data_set):
    """
    a class to hold the arrays used to map data to indexes of the nodes, edges
    or faces they are assigned to, if there is not that data on all of the objects
    do we ever need this?

    """
    ## fixme: do we want a special case for when there is data on all the nodes, edges, etc?

    def __init__(self, name, type='node', indexes=None, data=None):
        """
        create a data_set object

        name: the name of the data (depth, u_velocity, etc)

        type: the type of grid element: node, edge, or face the data is assigned to

        """
        self.name = name
        self.type = type # must be 'node', 'edge', of 'face' (eventually 'volume')
        if (indexes is None) ^  (data is None):
            raise ValueError("indexes and data both need to be either None or have values")
        if indexes is None:
            self.indexes = np.zeros((0,), dtype=IND_DT) 
        else:
            self.indexes = indexes
        if data is None:
            self.data = np.zeros((0,), dtype=np.float) # could be any data type
        else:
            self.data = data

    def check_consistent(self):
        """
        check if the indexes match the data, etc.
        """
        raise NotImplimentedError


class UGrid(object):
    """
    a basic class to hold an unstructred grid (triangular mesh)

    the internal structure mirrors the netcdf data standard.
    """

    def __init__(self,
                 nodes=None,
                 faces=None,
                 edges=None,
                 face_face_connectivity=None):
        """
        ugrid class -- holds, saves, etc an unstructured grid

        :param nodes=None : the coordinates of the nodes -- (NX2) float array
        :param faces=[] : the faces of the grid -- (NX3) integer array of indexes into the nodes array
        :param edges=[] : the edges of the grid -- (NX2) integer array of indexes into the nodes array

        often this is too much data to pass in a literal -- so usually
        specialized constructors will be used instead (load from file, etc.)
        """

        if nodes is None:
            self._nodes = np.zeros((0,2), dtype=NODE_DT)
        else:
            self._nodes = np.asarray(nodes, dtype=NODE_DT).reshape((-1, 2))

        if faces is None:
            self._faces = None
        else:
            self._faces = np.asarray(faces, dtype=IND_DT).reshape((-1, 3))

        if edges is None:
            self._edges = None
        else:
            self._edges = np.asarray(edges, dtype=IND_DT).reshape((-1, 2))

        if face_face_connectivity is None:
            self._face_face_connectivity = None
        else:
            self._face_face_connectivity = np.asarray(face_face_connectivity, dtype=IND_DT).reshape((-1, self.num_vertices))


        self._node_data = {}
        self._edge_data = {}
        self._face_data = {}

    def check_consistent(self):
        """
        check if the various data is consisent: the edges and faces reference
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
        self._nodes = np.asarray(nodes_coords, dtype=NODE_DT).reshape((-1, 2))
    @nodes.deleter
    def nodes(self):
        ## if there are no nodes, there can't be any faces or edges
        self._nodes = np.zeros((0,2), dtype=NODE_DT)
        self._edges = None
        self._faces = None

    @property
    def faces(self):
        return self._faces
    @faces.setter
    def faces(self, faces_indexes):
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
        self._faces = np.asarray(faces_indexes, dtype=IND_DT).reshape((-1, 3))
    @faces.deleter
    def faces(self):
        self._faces = None

    @property
    def edges(self):
        return self._edges
    @edges.setter
    def edges(self, edges_indexes):
        # room here to do consistency checking, etc.
        # for now -- simply make sure it's a numpy array
        self._edges = np.asarray(edges_indexes, dtype=IND_DT).reshape((-1, 2))
    @edges.deleter
    def edges(self):
        self._edges = None

    @property
    def face_face_connectivity(self):
        return self._face_face_connectivity
    @face_face_connectivity.setter
    def face_face_connectivity(self, face_face_connectivity):
        ## admore checking?
        face_face_connectivity = np.asarray(face_face_connectivity, dtype=IND_DT).reshape((-1, self.num_vertices))
    @face_face_connectivity.deleter
    def face_face_connectivity(self):
        self._face_face_connectivity = None



##fixme: repeated code here -- should these methods be combined?
    def set_node_data(self, name, data, indexes=None):
        if indexes is None:
            data = np.asarray(data)
            if not data.shape == (self.num_nodes,):
                raise ValueError("size of data should match number of nodes") # shape should match edges, data type can be anything
            self._node_data[name] = data
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            self._node_data[name][indexes] = data

    def get_node_data(self, name, indexes=None):
        if indexes is None:
            return self._node_data[name]
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            return self._node_data[name][indexes]

    def set_edge_data(self, name, data, indexes=None):
        if indexes is None:
            data = np.asarray(data)
            if not data.shape == (self.num_edges,):
                raise ValueError("size of data shold match number of edges") # shape should match edges, data type can be anything
            self._edge_data[name] = data
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            self._edge_data[name][indexes] = data

    def get_edge_data(self, name, indexes=None):
        if indexes is None:
            return self._edge_data[name]
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            return self._edge_data[name][indexes]

    def set_face_data(self, name, data, indexes=None):
        if indexes is None:
            data = np.asarray(data)
            if not data.shape == (self.num_faces,):
                raise ValueError("size of data shold match number of faces") # shape should match faces, data type can be anything
            self._face_data[name] = data
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            self._face_data[name][indexes] = data

    def get_face_data(self, name, indexes=None):
        if indexes is None:
            return self._face_data[name]
        else:
            indexes = np.array(indexes, dtype=IND_DT).reshape((-1,))
            return self._face_data[name][indexes]


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
        """        
        num_vertices = self.num_vertices
        num_faces = self.faces.shape[0]
        face_face = np.zeros( (num_faces, num_vertices), dtype=IND_DT )
        face_face += -1 # fill with -1

        # loop through all the triangles to find the matching edges:
        edges = {} # dict to store the edges in 
        for i, face in enumerate(self.faces):
            # loop through edges:
            for j in range(num_vertices):
                edge = (face[j-1], face[j])
                if edge[0] > edge[1]: # sort the node numbers
                    edge = (edge[1], edge[0]) 
                # see if it is already in there
                prev_edge = edges.pop(edge, None)
                if prev_edge is not None:
                    face_num, edge_num = prev_edge
                    face_face[i,j] = face_num
                    face_face[face_num, edge_num] = i
                else:
                    edges[edge] = (i, j)
        self._face_face_connectivity = face_face

    def build_edges(self):
        """
        builds the edges array: all the edges defined by the triangles

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


    def save_as_netcdf(self, filepath):
        """
        save the ugrid object as a netcdf file

        :param filepath: path to file you want o save to.
                         An existing one will be clobbered if it already exists.

        follows the convernsion established by the netcdf UGRID working group:

        http://publicwiki.deltares.nl/display/NETCDF/Deltares+CF+proposal+for+Unstructured+Grid+data+model

        """

        from netCDF4 import Dataset as ncDataset
        from netCDF4 import num2date, date2num
        # create a new netcdf file
        nclocal = ncDataset(filepath, mode="w", clobber=True)

        # dimensions:
        # nMesh2_node = 4 ; // nNodes
        # nMesh2_edge = 5 ; // nEdges
        # nMesh2_face = 2 ; // nFaces
        # nMesh2_face_links = 1 ; // nFacePairs


        nclocal.createDimension('num_nodes', len(self.nodes) )
        nclocal.createDimension('num_edges', len(self.edges) )
        nclocal.createDimension('num_faces', len(self.faces) )
        nclocal.createDimension('num_vertices', self.faces.shape[1] )
        nclocal.createDimension('two', 2)

        #mesh topology
        mesh = nclocal.createVariable('mesh', IND_DT, (), )
        mesh.cf_role = "mesh_topology" 
        mesh.long_name = "Topology data of 2D unstructured mesh" 
        mesh.topology_dimension = 2 
        mesh.node_coordinates = "node_lon node_lat" 
        mesh.face_node_connectivity = "mesh_face_nodes" 

        mesh.edge_node_connectivity = "mesh_edge_nodes"  ## attribute required if variables will be defined on edges
        mesh.edge_coordinates = "mesh_edge_lon mesh_edge_lat"  ## optional attribute (requires edge_node_connectivity)
        mesh.face_coordinates = "mesh_face_lon mesh_face_lat" ##  optional attribute
        mesh.face_edge_connectivity = "mesh_face_edges"  ## optional attribute (requires edge_node_connectivity)
        mesh.face_face_connectivity = "mesh_face_links"  ## optional attribute

        face_nodes = nclocal.createVariable("mesh_face_nodes",
                                            IND_DT,
                                            ('num_faces', 'num_vertices'),
                                            )
        face_nodes[:] = self.faces

        face_nodes.cf_role = "face_node_connectivity"
        face_nodes.long_name = "Maps every triangular face to its three corner nodes."
        face_nodes.start_index = 0 ;

        edge_nodes = nclocal.createVariable("mesh_edge_nodes",
                                            IND_DT,
                                            ('num_edges', 'two'),
                                            )
        edge_nodes[:] = self.edges

        edge_nodes.cf_role = "edge_node_connectivity"
        edge_nodes.long_name = "Maps every edge to the two nodes that it connects."
        edge_nodes.start_index = 0 ;

        node_lon = nclocal.createVariable('node_lon',
                                     self._nodes.dtype,
                                     ('num_nodes',),
                                     chunksizes=(len(self.nodes), ),
                                     #zlib=False,
                                     #complevel=0,
                                     )
        node_lon[:] = self.nodes[:,0]
        node_lon.standard_name = "longitude"
        node_lon.long_name = "Longitude of 2D mesh nodes."
        node_lon.units = "degrees_east"
        node_lat = nclocal.createVariable('node_lat',
                                     self._nodes.dtype,
                                     ('num_nodes',),
                                     chunksizes=(len(self.nodes), ),
                                     #zlib=False,
                                     #complevel=0,
                                     )
        node_lat[:] = self.nodes[:,1]
        node_lat.standard_name = "latitude"
        node_lat.long_name = "Latitude of 2D mesh nodes."
        node_lat.units = "degrees_north"

        # // Mesh node coordinates
        # double Mesh2_node_x(nMesh2_node) ;
        #         Mesh2_node_x:standard_name = "longitude" ;
        #         Mesh2_node_x:long_name = "Longitude of 2D mesh nodes." ;
        #         Mesh2_node_x:units = "degrees_east" ;
        nclocal.sync()







