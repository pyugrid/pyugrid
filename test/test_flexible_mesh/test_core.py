import os

import fiona
from numpy.ma import MaskedArray
import pytest
from shapely.geometry.base import BaseGeometry

from pyugrid.flexible_mesh.core import FlexibleMesh
from pyugrid.flexible_mesh.helpers import create_rtree_file, convert_multipart_to_singlepart, GeometryManager
from pyugrid.flexible_mesh.mpi import MPI_RANK, MPI_COMM
from test.test_flexible_mesh.base import AbstractFlexibleMeshTest


class TestFlexibleMesh(AbstractFlexibleMeshTest):

    def test_from_shapefile(self):
        shapefiles = {v: k for k, v in self.tdata_iter_shapefile_paths()}
        nfaces = {'disjoint': 2, 'single': 1, 'three': 3}
        out_nc_path = self.get_temporary_file_path('out.nc')
        out_shp_path = self.get_temporary_file_path('roundtrip.shp')

        keywords = dict(
            path=shapefiles.keys(),
            mesh_name=('_default', 'foobar'),
            # allow_multipolygons=(True, False),
            pack=(True, False),
            path_rtree=(None, '_create'),
        )

        for k in self.iter_product_keywords(keywords):

            name_uid = shapefiles[k.path]
            args = (k.path, name_uid)
            kwargs = {}
            if k.mesh_name != '_default':
                kwargs['mesh_name'] = k.mesh_name
            if k.path_rtree == '_create':
                rtree_path = self.get_temporary_file_path('rtree')
                rtree_gm = GeometryManager(name_uid, path=k.path)
                create_rtree_file(rtree_gm, rtree_path)
                kwargs['path_rtree'] = rtree_path
            kwargs['pack'] = k.pack

            res = FlexibleMesh.from_shapefile(*args, **kwargs)
            self.assertIsInstance(res, FlexibleMesh)

            # Test all faces are accounted for.
            found = False
            for kk, v in nfaces.items():
                if kk in k.path:
                    self.assertEqual(res.faces.shape[0], v)
                    found = True
                    break
            self.assertTrue(found)

            # Test writing the mesh back to netCDF.
            res.save_as_netcdf(out_nc_path)

            # Test reading the files back in.
            kwargs_from_ncfile = {'load_data': True}
            if k.mesh_name != '_default':
                kwargs_from_ncfile['mesh_name'] = k.mesh_name
            res2 = FlexibleMesh.from_ncfile(out_nc_path, **kwargs_from_ncfile)
            self.assertIsInstance(res2.faces, MaskedArray)
            self.assertEqual(res2.faces.shape, res.faces.shape)
            self.assertNumpyAll(res.faces, res2.faces, check_data=False, check_fill_value=False)

            # Test writing the read in file to shapefile.
            res2.save_as_shapefile(out_shp_path, face_uid_name=name_uid)
            self.assertShapefileGeometriesAlmostEqual(k.path, out_shp_path)
            uids = []
            for r in fiona.open(out_shp_path):
                uids.append(r['properties'][name_uid])
            self.assertTrue((res.data[name_uid].data == uids).all())

    def test_iter_records(self):
        path = self.get_temporary_file_path('out.shp')
        records, schema, name_uid = self.tdata_records_three
        self.write_fiona(path, records, schema)
        fm = FlexibleMesh.from_shapefile(path, name_uid)
        for shapely_only in [False, True]:
            for record in fm.iter_records(shapely_only=shapely_only):
                self.assertIn(name_uid, record['properties'])
                if shapely_only:
                    self.assertIsInstance(record['geom'], BaseGeometry)
                else:
                    self.assertIn('type', record['geometry'])

    @pytest.mark.mpi4py
    @pytest.mark.slow
    def test_save_as_netcdf_and_from_ncfile(self):
        path_nc = os.path.join(self.path_current_tmp, 'mesh2.nc')

        if MPI_RANK == 0:
            path_rtree = os.path.join(self.path_current_tmp, 'rtree')
            path_shp_single = os.path.join(self.path_current_tmp, 'shp_out.shp')
            path_shp_roundtrip = os.path.join(self.path_current_tmp, 'shp_out_roundtrip.shp')
            convert_multipart_to_singlepart(self.tdata_shapefile_path_state_boundaries, path_shp_single, start=5000)
            gm = GeometryManager('MID', path=path_shp_single)
            create_rtree_file(gm, path_rtree)
        else:
            path_shp_single = None
            path_rtree = None

        path_shp_single = MPI_COMM.bcast(path_shp_single, root=0)
        path_rtree = MPI_COMM.bcast(path_rtree, root=0)
        self.assertTrue(os.path.exists(path_rtree + '.idx'))

        fm = FlexibleMesh.from_shapefile(path_shp_single, 'MID', path_rtree=path_rtree)

        if MPI_RANK == 0:
            fm.save_as_shapefile(path_shp_roundtrip, face_uid_name='MID')
            self.assertShapefileGeometriesAlmostEqual(path_shp_single, path_shp_roundtrip)

            self.assertEqual(fm.data['MID'].data.shape, (134,))
            self.assertIsInstance(fm.faces, MaskedArray)
            self.assertEqual(fm.faces.shape[0], 134)
            self.assertEqual(fm.num_vertices, 1300)

            # Test converting to a netCDF, reading back in, and checking against the original.
            fm.save_as_netcdf(path_nc)
            fm2 = FlexibleMesh.from_ncfile(path_nc)
            self.assertNumpyAll(fm.faces, fm2.faces, check_data=False, check_fill_value=False)
