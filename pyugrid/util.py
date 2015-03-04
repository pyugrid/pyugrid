"""
miscellaneous util functions
"""
import numpy as np


epsilon = 1.e-5

def point_in_tri(face_points, point, return_weights=False):
    """ Calculates whether point is internal/external
    to element by comparing summed area of sub triangles with area of triangle element.
    """
    sub_tri_areas = np.zeros(3)
    sub_tri_areas[0] = _signed_area_tri(np.vstack((face_points[(0,1),:], point)))
    sub_tri_areas[1] = _signed_area_tri(np.vstack((face_points[(1,2),:], point)))
    sub_tri_areas[2] = _signed_area_tri(np.vstack((face_points[(0,2),:], point)))
    tri_area = _signed_area_tri(face_points)

    if abs(np.abs(sub_tri_areas).sum()-tri_area)/tri_area <= epsilon:
        if return_weights:
            raise NotImplementedError
            #weights = sub_tri_areas/tri_area
            #weights[1] = max(0., min(1., weights[1]))
            #weights[2] = max(0., min(1., weights[2]))
            #if (weights[0]+weights[1]>1):
            #    weights[2] = 0
            #    weights[1] = 1-weights[0]
            #else:
            #    weights[2] = 1-weights[0]-weights[1]
            #
            #return weights
        else:
            return True

    return False


def _signed_area_tri(points):
    """
    points : the coordinates of the triangle vertices -- (3X2) float array
    returns signed area of the triangle
    """

    x1, y1 = points[0]
    x2, y2 = points[1]
    x3, y3 = points[2]

    return(((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2)