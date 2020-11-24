from math import gcd as GCD

from tnorm.utilities.sage_types import Matrix
from tnorm.kernel.matrices import oriented_quads_mat
from tnorm.kernel.orientable import orient



def ends_embedded(qtons, TN_wrapper):
    pos_bdy, neg_bdy = TN_wrapper.boundary_slopes(qtons).values()
    for i in range(len(pos_bdy)):
        gcd_pos = gcd(pos_bdy[i][0],pos_bdy[i][1])
        gcd_neg = gcd(neg_bdy[i][0],neg_bdy[i][1])
        if gcd_pos != 0 and gcd_neg != 0:
            if not ( pos_bdy[i][0]/gcd_pos == - neg_bdy[i][0]/gcd_neg and pos_bdy[i][1]/gcd_pos == - neg_bdy[i][1]/gcd_neg ):
                return False
    return True

def is_embedded(qtons,TN_wrapper):
    W = TN_wrapper
    S = qtons
    T = W.triangulation()

    if not W.manifold_is_closed():
        pos_bdy, neg_bdy = W.boundary_slopes(S).values()
        if not ends_embedded(S,W):
            return False

    oriented_projection = orient(S)
    
    if oriented_projection != False:
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


def gcd(a,b):
    return abs(GCD(a,b))














