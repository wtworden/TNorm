#from sympy import Rational as QQ
from tnorm.sage_types import *

import tnorm.constants
from tnorm.kernel.matrices import oriented_quads_mat, intersection_mat, dot_product
from tnorm.kernel.regina_helpers import get_oriented_quads

SIGN_MATRIX = tnorm.constants.SIGN_MATRIX


### this maps an oriented spun normal surface to H_1(bdy M; \ZZ)
def qtons_to_H1bdy(oriented_spun_surface, TN_wrapper, quad_mat=None):
    s = oriented_spun_surface
    T = s.triangulation()

    M = TN_wrapper.manifold
    int_matrices = TN_wrapper._intersection_matrices

    sign_mat = SIGN_MATRIX  # SIGN_MATRIX encodes how the transverse orientation
    slopes = []                   # on a surface affects boundary orientations.
    for c in range(T.countCusps()):
        int_mat = int_matrices[c] # get the matrix that encodes meridian and longitude for cusp c.
        if quad_mat == None:
            quad_mat = oriented_quads_mat(s) # get the matrix that encodes quad coordinates for s.
        
        #iota_lambda = sum([Integer(int_mat[0][i]*sign_mat*quad_mat[i]) for i in range(T.size())])
        #iota_mu = sum([Integer(int_mat[1][i]*sign_mat*quad_mat[i]) for i in range(T.size())])
        iota_lambda=0
        iota_mu=0
        for i in range(T.size()):
            iota_lambda += Integer(int_mat[0][i]*sign_mat*quad_mat[i])
            iota_mu += Integer(int_mat[1][i]*sign_mat*quad_mat[i])

        ##### for sympy implementation
        #iota_lambda = QQ(sum([dot_product((int_mat[0].row(i)*sign_mat),quad_mat.row(i)) for i in range(T.size())]))
        #iota_mu = QQ(sum([dot_product((int_mat[1].row(i)*sign_mat),quad_mat.row(i)) for i in range(T.size())]))
        
        slopes.append((-iota_lambda,iota_mu))
    return slopes

### this maps an oriented spun normal surface to H_2(M,\partial M; \ZZ)
def homology_map(oriented_spun_surface, TN_wrapper=None):
    s = oriented_spun_surface
    slopes = qtons_to_H1bdy(s, TN_wrapper)
    coords = vector([slopes[i][1] for i in range(len(slopes))])
    return coords


#################--------------------------------------------------------------------



### this is the old version of Q_to_H1bdy()
#def intersection(cusp,oriented_spun_surface,M):
#    s = oriented_spun_surface
#    T = s.triangulation()
#    T.size()
#    mat = intersection_mat(M,T,cusp)
#    print(mat)
#    iota_lambda = 0
#    iota_mu = 0
#    for t in range(T.size()):
#        for v in range(4):
#            faces = [0,1,2,3]
#            faces.remove(v)
#            for f in faces:
#                verts = [0,1,2,3]
#                verts.remove(v)
#                verts.remove(f)
#                j,k = verts
#                iota_lambda += mat[0][t][4*v+f]*get_oriented_quads(s,t,j,k)*(-1)
#                iota_mu += mat[1][t][4*v+f]*get_oriented_quads(s,t,j,k)*(-1)
#                j,k = sorted([v,f])
#                iota_lambda += mat[0][t][4*v+f]*get_oriented_quads(s,t,j,k)
#                iota_mu += mat[1][t][4*v+f]*get_oriented_quads(s,t,j,k)
#    return (-iota_lambda,iota_mu)



