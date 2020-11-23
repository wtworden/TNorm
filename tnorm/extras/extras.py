from math import gcd as GCD


from tnorm.kernel import embedded, orientable, connected, euler
from tnorm.kernel.matrices import peripheral_curve_mats, oriented_quads_mat
from tnorm.kernel.boundary import signed_bdy_maps
import snappy
from tnorm.utilities.sage_types import *


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

def is_embedded(qtons):
    S = qtons
    T = S.triangulation()

    oriented_projection = orient(S)
    
    if oriented_projection != False:

        if not T.isClosed():
            pos_bdy, neg_bdy = oriented_boundary_slopes(S).values()
            if not ends_embedded(pos_bdy,neg_bdy):
                return False

        opp_oriented_projection = []
        for row in oriented_projection:
            opp_row = []
            for i in [0,2,4]:
                opp_row.append(row[i+1])
                opp_row.append(row[i])
            opp_oriented_projection.append(opp_row)

        q_mat = oriented_quads_mat(S)
        
        if oriented_projection == q_mat or Matrix(opp_oriented_projection) == q_mat:
            return True
    return False


def ends_embedded(pos_bdy, neg_bdy):
    for i in range(len(pos_bdy)):
        gcd_pos = gcd(pos_bdy[i][0],pos_bdy[i][1])
        gcd_neg = gcd(neg_bdy[i][0],neg_bdy[i][1])
        if gcd_pos != 0 and gcd_neg != 0:
            if not ( pos_bdy[i][0]/gcd_pos == - neg_bdy[i][0]/gcd_neg and pos_bdy[i][1]/gcd_pos == - neg_bdy[i][1]/gcd_neg ):
                return False
    return True


def gcd(a,b):
    return abs(GCD(a,b))