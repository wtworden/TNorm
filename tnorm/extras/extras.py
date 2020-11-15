
from tnorm.kernel import embedded, orientable, connected, euler
from tnorm.kernel.matrices import peripheral_curve_mats
from tnorm.kernel.boundary import signed_bdy_maps
import snappy


def orient(S):
    return orientable.orient(S)

def is_orientable(S):
    return orientable.is_orientable(S)

def euler_char(S):
    return euler.euler_char_(S,euler.solve_lin_gluing_eq(S.triangulation()))

def oriented_boundary_slopes(S):
    T = S.triangulation()
    M = snappy.Manifold(T.snapPea())
    pc_mats = peripheral_curve_mats(M,T)
    pos_bdy, neg_bdy = signed_bdy_maps(S, pc_mats, False, orient(S))
    return {'outward':pos_bdy, 'inward':neg_bdy}

def connected(S):
    pass