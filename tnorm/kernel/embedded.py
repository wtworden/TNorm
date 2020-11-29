from math import gcd as GCD
import itertools

from tnorm.utilities.sage_types import Matrix
from tnorm.kernel.matrices import oriented_quads_mat
from tnorm.kernel.orientable import orient



def ends_embedded(qtons, TN_wrapper):
    outward_bdy, inward_bdy = TN_wrapper.boundary_slopes(qtons).values()
    for i in range(len(outward_bdy)):
        gcd_outward = gcd(outward_bdy[i][0],outward_bdy[i][1])
        gcd_inward = gcd(inward_bdy[i][0],inward_bdy[i][1])
        if gcd_outward != 0 and gcd_inward != 0:
            if not ( outward_bdy[i][0]/gcd_outward == - inward_bdy[i][0]/gcd_inward and outward_bdy[i][1]/gcd_outward == - inward_bdy[i][1]/gcd_inward ):
                return False
    return True

def is_embedded(qtons,TN_wrapper):
    W = TN_wrapper
    S = qtons
    T = W.triangulation()

    if not W.manifold_is_closed():
        if not ends_embedded(S,W):
            return False

    oriented_projection = orient(S)

    if oriented_projection != False:
        k = len(oriented_projection)
        q_mat = oriented_quads_mat(S)

        orientations = list(itertools.product([1,-1], repeat=k))
        for nu in orientations:
            nu_oriented_proj = Matrix([[0 for i in range(6)] for j in range(T.size())])
            for i in range(k):
                if nu[i] == -1:
                    comp_orientation = opp_oriented_quads_mat(oriented_projection[i])
                else:
                    comp_orientation = oriented_projection[i]
                nu_oriented_proj += comp_orientation

            if nu_oriented_proj == q_mat:
                return True
  
    return False

def opp_oriented_quads_mat(oriented_quads_mat):
    opp_oriented_quads_mat = []
    for row in oriented_quads_mat:
        opp_row = []
        for i in [0,2,4]:
            opp_row.append(row[i+1])
            opp_row.append(row[i])
        opp_oriented_quads_mat.append(opp_row)
    return Matrix(opp_oriented_quads_mat)


def gcd(a,b):
    return abs(GCD(a,b))














