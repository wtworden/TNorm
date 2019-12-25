import snappy
import tnorm.constants
import tnorm.kernel.matrices
import tnorm.kernel.regina_helpers

from tnorm.kernel.matrices import oriented_quads_mat, intersection_mat
from tnorm.kernel.regina_helpers import regina_to_sage_int, regina_to_sage_mat

from tnorm.constants import ABS_SIGN_MATRIX

from tnorm.sage_types import *

#from tnorm.utilities import gcd
#from sympy import Rational as QQ



def get_quads(spun_surface,tet,v1,v2):
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_int(s.quads(tet,v2-int(1)))

## this map is almost the same as qtons_to_H1bdy(), except the SIGN_MATRIX is replaced by 
## SIGN_MATRIX.apply_map(abs). The result is the boundary slope for the embedded, possibly non-orientable,
## normal surface obtained by forgetting orientation. The reason this could be non-orientable is 
## that boundaries of the same slope that spin in oposite directions are cancelled. Thus the
## self-intersections coming from spinning are avoided, but at the cost of cutting and gluing
## in a way that does not respect orientation. This should return slopes which are (up to sign)
## the same as those returned by s.boundaryIntersections().
def bdy_slopes_unoriented_(oriented_spun_surface, TN_wrapper, quad_mat=None):
    s = oriented_spun_surface
    T = s.triangulation()

    M = TN_wrapper.manifold
    int_matrices = TN_wrapper._intersection_matrices

    sign_mat = ABS_SIGN_MATRIX # transverse orientation is ignored.
    slopes = [] 
    for c in range(T.countCusps()):
        int_mat = int_matrices[c] # get the matrix that encodes meridian and longitude for cusp c.
        if quad_mat == None:
            quad_mat = oriented_quads_mat(s) # get the matrix that encodes quad coordinates for s.
        iota_lambda = QQ(sum([int_mat[0][i]*sign_mat*quad_mat[i] for i in range(T.size())]))
        iota_mu = QQ(sum([int_mat[1][i]*sign_mat*quad_mat[i] for i in range(T.size())]))
        slopes.append((-iota_lambda,iota_mu))
    return slopes


def __boundary_slopes_(spun_surface, W=None):
    s = spun_surface
    T = s.triangulation()
    if W == None:
        M = snappy.Manifold(T.snapPea())
        int_matrices = [intersection_mat(M, T, cusp) for cusp in range(T.countCusps())]
    else:
        M = W.manifold
        int_matrices = W._intersection_matrices
    slopes = []
    for cusp in range(T.countCusps()):
        iota_lambda = 0
        iota_mu = 0
        mat = int_matrices[cusp]
        for t in range(T.size()):
            for v in range(4):
                faces = [0,1,2,3]
                faces.remove(v)
                for f in faces:
                    if 0 in [v,f]:
                        j,k = sorted([v,f])
                    else:
                        verts = [0,1,2,3]
                        verts.remove(v)
                        verts.remove(f)
                        j,k = verts
                    iota_lambda += mat[0][t][4*v+f]*get_quads(s,t,j,k)
                    iota_mu += mat[1][t][4*v+f]*get_quads(s,t,j,k)
        slopes.append((-iota_lambda,iota_mu))
    return slopes






