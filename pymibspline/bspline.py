import numpy as np
import itertools
import math

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


def findWeights(x, bins, spline_order, x_min = -1, x_max = -1):
    z = x2z(x, x_min, x_max)
    knots = knot_vector(bins, spline_order)
    ns = len(x)
    k = map(lambda e: basisFunction(e[0], spline_order, z[e[1]], knots, bins) , itertools.product(range(bins), range(ns)))
    return np.array(k).reshape(bins, ns)

def entropy1(w):
    (bins, ns) = w.shape
    h = w.sum(axis = 1)/ns
    return sum(map(lambda x: -x*math.log(x, 2) if x > 0 else 0  , h))

def entropy2(wx, wy):
    (garbage, ns) = wx.shape
    assert(ns == wy.shape[1])
    h = np.dot(wx, wy.T)/ns
    return sum(map(lambda x: -x*math.log(x, 2) if x > 0 else 0  , h.flatten()))

def mi(x, y, bins=5, so = 3, norm = True, negateMI = False):
    wx = findWeights(x, bins, so)
    wy = findWeights(y, bins, so)
    
    e1x = entropy1(wx)
    e1y = entropy1(wy)
    mi = e1x + e1y - entropy2(wx, wy)

    if norm:
        mix = 2*e1x - entropy2(wx, wx)
        miy = 2*e1y - entropy2(wy, wy)
        largerMI =  max(mix, miy)
        largerMI = 1 if largerMI == 0 else largerMI
        mi /= largerMI
    
    if negateMI and productMoment(x, y) < 0:
        mi = -mi
    
    return mi

def productMoment(x, y):
    return sum(map(lambda a: a[0]*a[1] ,zip(x, y))) * len(x) - sum(x) * sum(y)



