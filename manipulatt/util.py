# manipulatt - v2021.01.14
#Miscellaneous functions for manipulate.py
import numpy as np

try:
    from scipy.spatial import KDTree, Delaunay, distance_matrix
except ImportError:
    print('Please install scipy before proceeding')


def trim_filename(string): 
    return string.replace(" ", "") 
   
def _in_hull(p, hull):
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0
    
def centroid_coordinates(crds):
    return np.sum(crds[:, 0])/crds.shape[0], np.sum(crds[:, 1])/crds.shape[0], np.sum(crds[:, 2])/crds.shape[0]