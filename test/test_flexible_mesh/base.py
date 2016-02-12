from abc import ABCMeta
from collections import namedtuple, OrderedDict
from contextlib import contextmanager
from copy import deepcopy
import os
import shutil
import tempfile
from unittest import TestCase
import itertools

import netCDF4 as nc
import fiona
from shapely import wkt
from shapely.geometry import box, mapping, shape
import numpy as np


class AbstractFlexibleMeshTest(TestCase):
    __metaclass__ = ABCMeta
    key = 'ugrid_flexible_mesh'

    @property
    def tdata_records(self):
        mmo = OrderedDict([['disjoint', self.tdata_records_disjoint],
                           ['single', self.tdata_records_single],
                           ['three', self.tdata_records_three]])
        ret = OrderedDict()
        for k, v in mmo.items():
            ret[k] = {'records': v[0], 'schema': v[1], 'name_uid': v[2]}
        return ret

    @property
    def tdata_records_disjoint(self):
        record, schema, name_uid = self.tdata_records_single
        record.append({'geom': box(10, 20, 20, 30), 'properties': {'FID': 93}})
        return record, schema, name_uid

    @property
    def tdata_records_single(self):
        record = [{'geom': box(-1, -2, 1, 2), 'properties': {'FID': 91}}]
        schema = {'geometry': 'Polygon', 'properties': {'FID': 'int'}}
        return record, schema, 'FID'

    @property
    def tdata_records_three(self):
        square = 'Polygon ((-0.93061716886281376 0.32117641315631251, -0.39267504930065067 0.32571600910198478, -0.3994844432191591 -0.27351065572675382, -0.93288696683564987 -0.27124085775391771, -0.93061716886281376 0.32117641315631251))'
        triangle = 'Polygon ((0.01225455693106348 -0.26608971002281795, -0.3994844432191591 -0.27351065572675382, -0.39267504930065067 0.32571600910198478, 0.01225455693106348 -0.26608971002281795))'
        ngon = 'Polygon ((-0.20422637505866323 -0.26999143718892821, 0.01225455693106348 -0.26608971002281795, 0.27769939343655725 -0.30936875945306097, 0.27769939343655725 -0.57481359595855341, -0.07718881189143811 -0.64694534500895906, -0.26761662938450881 -0.47382914728798564, -0.20422637505866323 -0.26999143718892821))'
        records = [{'geom': wkt.loads(square), 'properties': {'SPECIAL': 100}},
                   {'geom': wkt.loads(triangle), 'properties': {'SPECIAL': 101}},
                   {'geom': wkt.loads(ngon), 'properties': {'SPECIAL': 102}}]
        return records, {'geometry': 'Polygon', 'properties': {'SPECIAL': 'int'}}, 'SPECIAL'

    @property
    def tdata_shapefile_path_state_boundaries(self):
        return os.path.join(self.path_files, 'state_boundaries', 'state_boundaries.shp')

    @property
    def path_files(self):
        import test
        path = os.path.split(test.__file__)[0]
        path = os.path.join(path, 'files')
        return path

    def assertNumpyAll(self, arr1, arr2, check_fill_value=True, check_arr_dtype=True, check_arr_type=True,
                       check_data=True):
        """
        Asserts arrays are equal according to the test criteria.

        :param arr1: An array to compare.
        :type arr1: :class:`numpy.ndarray`
        :param arr2: An array to compare.
        :type arr2: :class:`numpy.ndarray`
        :param bool check_fill_value_dtype: If ``True``, check if fill values are equal.
        :param bool check_arr_dtype: If ``True``, check the data types of the arrays are equal.
        :param bool check_arr_type: If ``True``, check the types of the incoming arrays:

        >>> type(arr1) == type(arr2)

        :param check_data: If ``True``, check that ``data`` attributes of masked arrays are equal.
        :type check_data: bool
        :raises: AssertionError
        """

        if check_arr_type:
            self.assertEqual(type(arr1), type(arr2))
        self.assertEqual(arr1.shape, arr2.shape)
        if check_arr_dtype:
            self.assertEqual(arr1.dtype, arr2.dtype)
        if isinstance(arr1, np.ma.MaskedArray) or isinstance(arr2, np.ma.MaskedArray):
            if check_data:
                self.assertTrue(np.all(arr1.data == arr2.data))
            else:
                self.assertTrue(np.all(arr1 == arr2))
            self.assertTrue(np.all(arr1.mask == arr2.mask))
            if check_fill_value:
                self.assertEqual(arr1.fill_value, arr2.fill_value)
        else:
            self.assertTrue(np.all(arr1 == arr2))

    def assertPolygonAlmostEqual(self, a, b):
        self.assertEqual(type(a), type(b))
        self.assertEqual(a.bounds, b.bounds)
        self.assertAlmostEqual(a.area, b.area)
        a_coords = np.array(a.exterior.coords)
        b_coords = np.array(b.exterior.coords)
        self.assertEqual(len(a_coords), len(b_coords))
        self.assertAlmostEqual(a_coords.mean(), b_coords.mean())

    def assertShapefileGeometriesAlmostEqual(self, lhs, rhs):
        with fiona.open(lhs) as source:
            polygons_original = [shape(e['geometry']) for e in source]

        with fiona.open(rhs) as s1:
            self.assertEqual(len(list(s1)), len(polygons_original))
            for r1 in s1:
                found = False
                g1 = shape(r1['geometry'])
                for g2 in polygons_original:
                    try:
                        self.assertPolygonAlmostEqual(g1, g2)
                        found = True
                        break
                    except AssertionError:
                        continue
                self.assertTrue(found)

    def get_temporary_file_path(self, fn):
        return os.path.join(self.path_current_tmp, fn)

    @staticmethod
    @contextmanager
    def nc_scope(path, mode='r', format=None):
        """
        Provide a transactional scope around a :class:`netCDF4.Dataset` object.

        >>> with nc_scope('/my/file.nc') as ds:
        >>>     print ds.variables

        :param str path: The full path to the netCDF dataset.
        :param str mode: The file mode to use when opening the dataset.
        :param str format: The NetCDF format.
        :returns: An open dataset object that will be closed after leaving the ``with`` statement.
        :rtype: :class:`netCDF4.Dataset`
        """

        kwds = {'mode': mode}
        if format is not None:
            kwds['format'] = format

        ds = nc.Dataset(path, **kwds)
        try:
            yield ds
        finally:
            ds.close()

    def iter_product_keywords(self, keywords, as_namedtuple=True):
        return itr_products_keywords(keywords, as_namedtuple=as_namedtuple)

    def tdata_iter_shapefile_paths(self):
        for k, v in self.tdata_records.items():
            path = self.get_temporary_file_path(k + '.shp')
            self.write_fiona(path, v['records'], v['schema'])
            yield v['name_uid'], path

    def write_fiona(self, path, records, schema, driver='ESRI Shapefile'):
        with fiona.open(path, 'w', schema=schema, driver=driver) as sink:
            for r in records:
                r['geometry'] = mapping(r['geom'])
                sink.write(r)

    def setUp(self):
        self.path_current_tmp = tempfile.mkdtemp(prefix='{0}_test_'.format(self.key))

    def shortDescription(self):
        return None

    def tearDown(self):
        shutil.rmtree(self.path_current_tmp)


def attr(*args, **kwargs):
    """
    Decorator that adds attributes to classes or functions for use with the Attribute (-a) plugin.

    http://nose.readthedocs.org/en/latest/plugins/attrib.html
    """

    def wrap_ob(ob):
        for name in args:
            setattr(ob, name, True)
        for name, value in kwargs.items():
            setattr(ob, name, value)
        return ob

    return wrap_ob


def itr_row(key, sequence):
    for element in sequence:
        yield ({key: element})


def itr_products_keywords(keywords, as_namedtuple=False):
    if as_namedtuple:
        yld_tuple = namedtuple('ITesterKeywords', keywords.keys())

    iterators = [itr_row(ki, vi) for ki, vi in keywords.items()]
    for dictionaries in itertools.product(*iterators):
        yld = {}
        for dictionary in dictionaries:
            yld.update(dictionary)
        if as_namedtuple:
            yld = yld_tuple(**yld)
        yield deepcopy(yld)
