import sys
import numpy as np
import itertools
import math
import collections
import _c_bsplinemi
from LabeledMat import LabeledMat

def knot_vector(bins, spline_order):
    internal_points = bins - spline_order + 1
    v = [[0]*spline_order , range(1, internal_points)  , [internal_points]*spline_order]
    return [x/float(internal_points) for vv in v for x in vv]

def x2z(x, x_min = -1, x_max = -1):
    if x_min == -1 and x_max == -1:
        x_min = min(x)
        x_max = max(x)
    if x_min == x_max:
        x_max = x_min + 1
    return map(lambda xx: (xx - x_min) / float(x_max - x_min)  , x)

def basisFunction(i, p, t, kVector, bins):
    if p == 1:
        if (t >= kVector[i] and t < kVector[i+1] and kVector[i] < kVector[i+1]) or (abs(t - kVector[i+1]) < 1E-10 and i+1 == bins):
            return 1
        else:
            return 0
    d1 = kVector[i+p-1] - kVector[i]
    n1 = t - kVector[i]
    d2 = kVector[i+p] - kVector[i+1]
    n2 = kVector[i+p] - t
    if d1 < 1E-10 and d2 < 1E-10:
        return 0
    elif d1 < 1E-10:
        e1 = 0
        e2 = n2 / d2 * basisFunction(i+1, p-1, t, kVector, bins)
    elif d2 < 1E-10:
        e2 = 0
        e1 = n1 / d1 * basisFunction(i, p-1, t, kVector, bins)
    else:
        e1 = n1 / d1 * basisFunction(i, p-1, t, kVector, bins)
        e2 = n2 / d2 * basisFunction(i+1, p-1, t, kVector, bins)
    if e1 + e2 < 0:
        return 0
    return e1+e2 

def find_weights(x, bins = 6, so = 3, x_min = -1, x_max = -1):
    if not isinstance(x, collections.Iterable):
        print >> sys.stderr, "ERROR: input vector is not iterable!"
        raise
    n = len(x)
    w = _c_bsplinemi.find_weights(x, bins, spline_order)
    return w.reshape((bins, n))

def entropy(x, bins=6, so=3):
    if not isinstance(x, collections.Iterable):
        print >> sys.stderr, "ERROR: input vector is not iterable!"
        raise
    return _c_bsplinemi.entropy(x, bins, so)

def joint_entropy(x, y, bins=6, so=3):
    if not isinstance(x, collections.Iterable) or not isinstance(y, collections.Iterable):
        print >> sys.stderr, "ERROR: input vectors are not both iterable!"
        raise

    if len(x) != len(y):
        print >> sys.stderr, "ERROR: two vectors must be of same length!"
        raise

    return _c_bsplinemi.joint_entropy(x, y, bins, so)

def mi(x, y, bins=6, so = 3, norm = True, negateMI = False):
    if not isinstance(x, collections.Iterable) or not isinstance(y, collections.Iterable):
        print >> sys.stderr, "ERROR: input vectors are not both iterable!"
        raise

    if len(x) != len(y):
        print >> sys.stderr, "ERROR: two vectors must be of same length!"
        raise

    return _c_bsplinemi.mi(x, y, bins, so, norm, negateMI)

def all_mi(X, vec, bins=6, so = 3, norm=True, negateMI=True):
    if X.__class__.__name__ != 'LabeledMat':
        print >> sys.stderr, "ERROR: input matrix must be LabeledMat!"
        raise

    if not isinstance(vec, collections.Iterable):
        print >> sys.stderr, "ERROR: input vector must be iterable!"
        raise

    if X.ncol != len(vec):
        print >> sys.stderr, "ERROR: two vectors must be of same length!"
        raise

    mis = _c_bsplinemi.all_mi(X.data, vec, bins, so, norm, negateMI)
    return dict(zip(X.rownames, mis))



