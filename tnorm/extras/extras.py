from math import gcd as GCD
import itertools

from tnorm.kernel import embedded, orientable, connected, euler
from tnorm.kernel.matrices import peripheral_curve_mats, oriented_quads_mat, quads_mat
from tnorm.kernel.boundary import inward_oriented_bdy, outward_oriented_bdy, unoriented_spinning_slopes
import snappy

from tnorm.utilities.sage_types import Matrix

from tnorm.utilities.regina_helpers import regina_to_sage_int


def orient(S):
    return orientable.orient(S)

def is_orientable(S):
    return orientable.is_orientable(S)

def euler_char(S):
    T = S.triangulation()
    assert T.isIdeal(), 'the triangulation must be ideal.'
    return euler.euler_char_(S,euler.solve_lin_gluing_eq(T))

def spinning_slopes(S):
    q_mat = quads_mat(S)
    T = S.triangulation()
    M = snappy.Manifold(T.snapPea())
    pc_mats = peripheral_curve_mats(M,T)
    slopes = unoriented_spinning_slopes(S, pc_mats, q_mat)
    return slopes

def num_boundary_comps(S):
    slopes = spinning_slopes(S)
    return sum([gcd(slope[0],slope[1]) for slope in slopes])

def oriented_boundary_slopes(S, orientation=None):
    if orientation == None:
        oriented_quads_mat = sum(orient(S))
    else:
        oriented_quads_mat = orientation
    assert oriented_quads_mat != False, 'the surface is not orientable.'
    T = S.triangulation()
    M = snappy.Manifold(T.snapPea())
    pc_mats = peripheral_curve_mats(M,T)
    out_bdy = outward_oriented_bdy(S, pc_mats, oriented_quads_mat, False)
    in_bdy = inward_oriented_bdy(S, pc_mats, oriented_quads_mat, False)

    return {'outward':out_bdy, 'inward':in_bdy}, oriented_quads_mat

def is_connected(S):
    return connected.is_connected(S)

def connected_components(S):
    return connected.connected_components(S)

def haken_sum(R,S):
    if R.locallyCompatible(S):
        return connected.sum(R,S)
    else:
        return "these surfaces are not locally compatible!"


def gcd(a,b):
    return abs(GCD(a,b))