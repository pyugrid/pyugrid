from collections import deque, OrderedDict
from copy import copy
import os

import numpy as np
import fiona
from shapely.geometry import shape, mapping, Polygon
from shapely.geometry.base import BaseMultipartGeometry
from shapely.geometry.polygon import orient

from pyugrid.flexible_mesh.constants import PYUGRID_LINK_ATTRIBUTE_NAME
from pyugrid.flexible_mesh.geom_cabinet import GeomCabinetIterator
from pyugrid.flexible_mesh.mpi import MPI_RANK, create_slices, MPI_COMM, hgather, vgather, MPI_SIZE
from pyugrid.flexible_mesh.logging_ugrid import log
from pyugrid.flexible_mesh.spatial_index import SpatialIndex


def convert_multipart_to_singlepart(path_in, path_out, new_uid_name=PYUGRID_LINK_ATTRIBUTE_NAME, start=0):
    """
    Convert a vector GIS file from multipart to singlepart geometries. The function copies all attributes and
    maintains the coordinate system.

    :param str path_in: Path to the input file containing multipart geometries.
    :param str path_out: Path to the output file.
    :param str new_uid_name: Use this name as the default for the new unique identifier.
    :param int start: Start value for the new unique identifier.
    """

    with fiona.open(path_in) as source:
        source.meta['schema']['properties'][new_uid_name] = 'int'
        with fiona.open(path_out, mode='w', **source.meta) as sink:
            for record in source:
                geom = shape(record['geometry'])
                if isinstance(geom, BaseMultipartGeometry):
                    for element in geom:
                        record['properties'][new_uid_name] = start
                        record['geometry'] = mapping(element)
                        sink.write(record)
                        start += 1
                else:
                    record['properties'][new_uid_name] = start
                    sink.write(record)
                    start += 1


def iter_edge_nodes(idx_nodes):
    for ii in range(len(idx_nodes)):
        try:
            yld = (idx_nodes[ii], idx_nodes[ii + 1])
        # the last node pair requires linking back to the first node
        except IndexError:
            yld = (idx_nodes[-1], idx_nodes[0])
        yield yld


def get_variables(gm, pack=False):
    """
    :param gm: The geometry manager containing geometries to convert to mesh variables.
    :type gm: :class:`pyugrid.flexible_mesh.helpers.GeometryManager`
    :param pack: If ``True``, de-deduplicate shared coordinates.
    :type pack: bool
    :returns: A tuple of arrays with index locations corresponding to:

    ===== ================ =============================
    Index Name             Type
    ===== ================ =============================
    0     face_nodes       :class:`numpy.ma.MaskedArray`
    1     face_edges       :class:`numpy.ma.MaskedArray`
    2     edge_nodes       :class:`numpy.ndarray`
    3     node_x           :class:`numpy.ndarray`
    4     node_y           :class:`numpy.ndarray`
    5     face_links       :class:`numpy.ndarray`
    6     face_ids         :class:`numpy.ndarray`
    7     face_coordinates :class:`numpy.ndarray`
    ===== ================ =============================

    Information on individual variables may be found here: https://github.com/ugrid-conventions/ugrid-conventions/blob/9b6540405b940f0a9299af9dfb5e7c04b5074bf7/ugrid-conventions.md#2d-flexible-mesh-mixed-triangles-quadrilaterals-etc-topology

    :rtype: tuple (see table for array types)
    :raises: ValueError
    """

    if len(gm) < MPI_SIZE:
        raise ValueError('the number of geometries must be greater than or equal to the number of processes')

    result = get_face_variables(gm)

    if MPI_RANK == 0:
        face_links, nmax_face_nodes, face_ids, face_coordinates = result
    else:
        return

    # the number of faces
    n_face = len(gm)
    face_nodes = np.ma.array(np.zeros((n_face, nmax_face_nodes), dtype=np.int32), mask=True)
    # the edge mapping has the same shape as the node mapping
    face_edges = np.zeros_like(face_nodes)

    # holds the start and end nodes for each edge
    edge_nodes = deque()

    # holds x and y node values
    node_x = np.array([])
    node_y = np.array([])
    node_indices = np.array([])

    # flag to indicate if this is the first face encountered
    first = True
    # global point index counter
    idx_point = 0
    # global edge index counter
    idx_edge = 0
    # loop through each polygon
    for idx_face_nodes, (fid, record) in enumerate(gm.iter_records(return_uid=True)):
        coordinates = np.array(record['geom'].exterior.coords)
        # Assert last coordinate is repeated for each polygon.
        assert coordinates[0].tolist() == coordinates[-1].tolist()
        # Remove this repeated coordinate.
        coordinates = coordinates[0:-1, :]

        # just load everything if this is the first polygon
        if first:
            # store the point values.
            for ii in range(coordinates.shape[0]):
                coordinates_row = coordinates[ii, :]
                # store the x and y coordinates
                node_x = np.hstack((node_x, [coordinates_row[0]]))
                node_y = np.hstack((node_y, [coordinates_row[1]]))
                # increment the point index
                idx_point += 1
            # map the node identifiers for the face
            range_point_idx = range(0, idx_point)
            face_nodes[idx_face_nodes, 0:idx_point] = range_point_idx
            node_indices = np.hstack((node_indices, range_point_idx))
            # construct the edges. compress the node slice to remove any masked values at the tail.
            for start_node_idx, end_node_idx in iter_edge_nodes(face_nodes[idx_face_nodes, :].compressed()):
                edge_nodes.append((start_node_idx, end_node_idx))
                idx_edge += 1
            # map the edges to faces
            face_edges[idx_face_nodes, 0:idx_edge] = range(0, idx_edge)
            # switch the loop flag to indicate the first face has been dealt with
            first = False
        else:
            # holds new node coordinates for the face
            new_face_nodes = deque()

            # Only search neighboring faces if there are face links.
            if face_links is not None:
                neighbor_face_indices = face_links == idx_face_nodes
                neighbor_face_indices = np.logical_or(neighbor_face_indices[:, 0], neighbor_face_indices[:, 1])
                neighbor_face_indices = face_links[neighbor_face_indices, :]
                neighbor_face_indices = neighbor_face_indices.flatten()[
                    np.where(neighbor_face_indices.flatten() != idx_face_nodes)[0]]
                neighbor_face_indices = np.unique(neighbor_face_indices)
            else:
                neighbor_face_indices = np.array([])

            node_indices, node_x, node_y, idx_point = get_mapped_xy(face_nodes, node_indices, node_x, node_y,
                                                                    coordinates, idx_point, neighbor_face_indices,
                                                                    new_face_nodes, pack=pack)
            # map the node indices for the face
            face_nodes[idx_face_nodes, 0:len(new_face_nodes)] = new_face_nodes
            # find and map the edges
            idx_edge = get_mapped_edges(edge_nodes, face_edges, face_nodes, idx_face_nodes, idx_edge, pack=pack)

    return face_nodes, face_edges, np.array(edge_nodes, dtype=np.int32), node_x, node_y, face_links, face_ids, \
           face_coordinates


def get_mapped_xy(face_nodes, node_indices, node_x, node_y, coordinates, idx_point, neighbor_face_indices,
                  new_face_nodes, pack=False):
    for ii in range(coordinates.shape[0]):
        # logic flag to indicate if the point has been found
        found = False
        coordinates_row = coordinates[ii, :]
        if pack:
            # search the neighboring faces for matching nodes
            for neighbor_face_index in neighbor_face_indices.flat:
                # search over the neighboring face's nodes
                neighbor_face_node_indices = face_nodes[neighbor_face_index, :].compressed()
                # for neighbor_face_node_index in compressed:
                x_equal = np.isclose(coordinates_row[0], node_x[neighbor_face_node_indices])
                y_equal = np.isclose(coordinates_row[1], node_y[neighbor_face_node_indices])
                is_equal = np.logical_and(x_equal, y_equal)
                if np.any(is_equal):
                    new_face_nodes.extend(node_indices[neighbor_face_node_indices][is_equal].tolist())
                    # point is found, no need to continue with loop
                    found = True
                    break
        # add the new node if it has not been found
        if not found:
            # add the coordinates of the new point
            node_x = np.hstack((node_x, coordinates_row[0]))
            node_y = np.hstack((node_y, coordinates_row[1]))
            # append the index of this new point
            new_face_nodes.append(idx_point)
            # increment the point index
            node_indices = np.hstack((node_indices, idx_point))
            idx_point += 1

    return node_indices, node_x, node_y, idx_point


def get_mapped_edges(edge_nodes, face_edges, face_nodes, idx_face_nodes, idx_edge, pack=False):
    new_face_edges = deque()
    for start_node_idx, end_node_idx in iter_edge_nodes(face_nodes[idx_face_nodes, :].compressed()):
        # flag to indicate if edge has been found
        found_edge = False
        # search existing edge-node combinations accounting for ordering
        if pack:
            found_edge = get_found_edge(edge_nodes, end_node_idx, found_edge, new_face_edges, start_node_idx)
        if not found_edge:
            edge_nodes.append((start_node_idx, end_node_idx))
            new_face_edges.append(idx_edge)
            idx_edge += 1
    # update the face-edge mapping
    face_edges[idx_face_nodes, 0:len(new_face_edges)] = new_face_edges
    return idx_edge


def get_found_edge(edge_nodes, end_node_idx, found_edge, new_face_edges, start_node_idx):
    for idx_edge_nodes, edge_nodes in enumerate(edge_nodes):
        # swap the node ordering
        if edge_nodes == (start_node_idx, end_node_idx) or edge_nodes == (end_node_idx, start_node_idx):
            new_face_edges.append(idx_edge_nodes)
            found_edge = True
            break
    return found_edge


def get_found_edge_indexes(edge_nodes, end_node_idx, start_node_idx):
    sums = np.array(edge_nodes).sum(axis=1)
    close = np.isclose(sums, np.sum([start_node_idx, end_node_idx]))
    indexes = np.where(close)[0]
    return indexes


def iter_touching(si, gm, shapely_object):
    select_uid = list(si.iter_rtree_intersection(shapely_object))
    select_uid.sort()
    for uid_target, record_target in gm.iter_records(return_uid=True, select_uid=select_uid):
        if shapely_object.touches(record_target['geom']):
            yield uid_target


def get_face_variables(gm):

    if MPI_RANK == 0:
        sections = create_slices(len(gm))
    else:
        sections = None

    section = MPI_COMM.scatter(sections, root=0)

    # Create a spatial index to find touching faces.
    si = gm.get_spatial_index()

    Mesh2_ids = np.zeros(section[1] - section[0], dtype=np.int32)
    assert Mesh2_ids.shape[0] > 0

    Mesh2_face_links = deque()
    nMaxMesh2_face_nodes = 0
    face_coordinates = deque()

    log.debug('section={0} (rank={1})'.format(section, MPI_RANK))

    for ctr, (uid_source, record_source) in enumerate(gm.iter_records(return_uid=True, slice=section)):
        Mesh2_ids[ctr] = uid_source
        ref_object = record_source['geom']

        # Get representative points for each polygon.
        face_coordinates.append(np.array(ref_object.representative_point()))

        # For polygon geometries the first coordinate is repeated at the end of the sequence. UGRID clients do not want
        # repeated coordinates (i.e. ESMF).
        ncoords = len(ref_object.exterior.coords) - 1
        if ncoords > nMaxMesh2_face_nodes:
            nMaxMesh2_face_nodes = ncoords

        for uid_target in iter_touching(si, gm, ref_object):
            # If the objects only touch they are neighbors and may share nodes.
            Mesh2_face_links.append([uid_source, uid_target])

    # convert to numpy array for faster comparisons
    Mesh2_face_links = np.array(Mesh2_face_links, dtype=np.int32)

    Mesh2_ids = MPI_COMM.gather(Mesh2_ids, root=0)
    nMaxMesh2_face_nodes = MPI_COMM.gather(nMaxMesh2_face_nodes)
    Mesh2_face_links = MPI_COMM.gather(Mesh2_face_links)
    face_coordinates = MPI_COMM.gather(np.array(face_coordinates))

    if MPI_RANK == 0:
        Mesh2_ids = hgather(Mesh2_ids)
        face_coordinates = vgather(face_coordinates)

        nMaxMesh2_face_nodes = max(nMaxMesh2_face_nodes)

        # For single element conversion, the links should be empty (i.e. no neighbors).
        try:
            Mesh2_face_links = vgather(Mesh2_face_links)
        except IndexError:
            if Mesh2_face_links[0].shape[0] == 0:
                Mesh2_face_links = None
            else:
                raise
        else:
            new_face_links = Mesh2_face_links.copy()
            # Convert face unique identifiers to indices.
            for ctr, mid in enumerate(Mesh2_ids.flat):
                for column_index in [0, 1]:
                    w = np.where(Mesh2_face_links[:, column_index] == mid)[0]
                    new_face_links[w, column_index] = ctr
            Mesh2_face_links = new_face_links

        return Mesh2_face_links, nMaxMesh2_face_nodes, Mesh2_ids, face_coordinates


def flexible_mesh_to_fiona(out_path, face_nodes, node_x, node_y, crs=None, driver='ESRI Shapefile',
                           indices_to_load=None, face_uid=None):

    if face_uid is None:
        properties = {}
    else:
        properties = {face_uid.name: 'int'}

    schema = {'geometry': 'Polygon', 'properties': properties}
    with fiona.open(out_path, 'w', driver=driver, crs=crs, schema=schema) as f:
        for feature in iter_records(face_nodes, node_x, node_y, indices_to_load=indices_to_load, datasets=[face_uid]):
            feature['properties'][face_uid.name] = int(feature['properties'][face_uid.name])
            f.write(feature)
    return out_path


def iter_records(face_nodes, node_x, node_y, indices_to_load=None, datasets=None, shapely_only=False):
    if indices_to_load is None:
        feature_indices = range(face_nodes.shape[0])
    else:
        feature_indices = indices_to_load

    for feature_idx in feature_indices:
        coordinates = deque()

        current_face_node = face_nodes[feature_idx, :]
        try:
            nodes = current_face_node.compressed()
        except AttributeError:
            # Likely not a masked array.
            nodes = current_face_node.flatten()

        # Construct the geometry object by collecting node coordinates using indicies stored in "nodes".
        for node_idx in nodes:
            coordinates.append((node_x[node_idx], node_y[node_idx]))
        polygon = Polygon(coordinates)

        # Collect properties if datasets are passed.
        properties = OrderedDict()
        if datasets is not None:
            for ds in datasets:
                properties[ds.name] = ds.data[feature_idx]
        feature = {'id': feature_idx, 'properties': properties}

        # Add coordinates or shapely objects depending on parameters.
        if shapely_only:
            feature['geom'] = polygon
        else:
            feature['geometry'] = mapping(polygon)

        yield feature


def create_rtree_file(gm, path):
    """
    :param gm: Target geometries to index.
    :type gm: :class:`pyugrid.flexible_mesh.helpers.GeometryManager`
    :param path: Output path for the serialized spatial index. See http://toblerity.org/rtree/tutorial.html#serializing-your-index-to-a-file.
    """

    si = SpatialIndex(path=path)
    for uid, record in gm.iter_records(return_uid=True):
        si.add(uid, record['geom'])


class GeometryManager(object):
    """
    Provides iteration, validation, and other management routines for collecting vector geometries from record lists or
    flat files.
    """

    def __init__(self, name_uid, path=None, records=None, path_rtree=None):
        if path_rtree is not None:
            assert os.path.exists(path_rtree + '.idx')

        self.path = path
        self.path_rtree = path_rtree
        self.name_uid = name_uid
        self.records = copy(records)

        self._has_provided_records = False if records is None else True

    def __len__(self):
        if self.records is None:
            ret = len(GeomCabinetIterator(path=self.path))
        else:
            ret = len(self.records)
        return ret

    def get_spatial_index(self):
        si = SpatialIndex(path=self.path_rtree)
        # Only add new records to the index if we are working in-memory.
        if self.path_rtree is None:
            for uid, record in self.iter_records(return_uid=True):
                si.add(uid, record['geom'])
        return si

    def iter_records(self, return_uid=False, select_uid=None, slice=None):
        # Use records attached to the object or load records from source data.
        to_iter = self.records or self._get_records_(select_uid=select_uid, slice=slice)

        if self.records is not None and slice is not None:
            to_iter = to_iter[slice[0]:slice[1]]

        for ctr, record in enumerate(to_iter):
            if self._has_provided_records and 'geom' not in record:
                record['geom'] = shape(record['geometry'])
                # Only use the geometry objects from here. Maintaining the list of coordinates is superfluous.
                record.pop('geometry')
            self._validate_update_record_(record)

            # Counter-clockwise orientations required by clients such as ESMF Mesh regridding.
            obj = record['geom']
            if not obj.exterior.is_ccw:
                record['geom'] = orient(obj)

            if return_uid:
                uid = record['properties'][self.name_uid]
                yld = (uid, record)
            else:
                yld = record
            yield yld

    def _get_records_(self, select_uid=None, slice=slice):
        gi = GeomCabinetIterator(path=self.path, uid=self.name_uid, select_uid=select_uid, slice=slice)
        return gi

    def _validate_update_record_(self, record):
        geom = record['geom']

        # This should happen before any buffering. The buffering check may result in a single polygon object.
        if isinstance(geom, BaseMultipartGeometry):
            msg = 'Only singlepart geometries allowed. Perhaps "ugrid.convert_multipart_to_singlepart" would be useful?'
            raise ValueError(msg)

        try:
            assert geom.is_valid
        except AssertionError:
            geom = geom.buffer(0)
            assert geom.is_valid
            record['geom'] = geom