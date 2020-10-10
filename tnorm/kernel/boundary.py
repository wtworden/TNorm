import snappy
import tnorm.constants
import tnorm.kernel.matrices
import tnorm.kernel.regina_helpers

from tnorm.kernel.matrices import oriented_quads_mat, peripheral_curve_mats, quad_bdy_mat, intersection_mat, pos_intersection_mat, neg_intersection_mat  
from tnorm.kernel.regina_helpers import regina_to_sage_int, regina_to_sage_mat


from tnorm.sage_types import *

#from tnorm.utilities import gcd
#from sympy import Rational as QQ



def get_quads(spun_surface,tet,v1,v2):
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_int(s.quads(tet,v2-int(1)))


### this maps an oriented spun normal surface to H_1(bdy M; \ZZ)
def signed_bdy_maps(oriented_spun_surface, TN_wrapper, spinning=False):
    s = oriented_spun_surface
    T = s.triangulation()

    M = TN_wrapper.manifold()
    periph_mats = TN_wrapper._peripheral_curve_mats

    pos_intx_mat = pos_intersection_mat()
    neg_intx_mat =  neg_intersection_mat()
    if spinning:
        neg_intx_mat = neg_intx_mat.apply_map(abs)
    pos_slopes = []                   # on a surface affects boundary orientations.
    neg_slopes = []
    quads_mat = oriented_quads_mat(s)
    for c in range(T.countCusps()):
        periph_mat = periph_mats[c] # get the matrix that encodes meridian and longitude for cusp c.
        
        iota_lambda_pos=0
        iota_mu_pos=0
        iota_lambda_neg=0
        iota_mu_neg=0
        for i in range(T.size()):
            iota_lambda_pos += Integer(periph_mat[0][i]*pos_intx_mat*quads_mat[i])
            iota_mu_pos += Integer(periph_mat[1][i]*pos_intx_mat*quads_mat[i])
            iota_lambda_neg += Integer(periph_mat[0][i]*neg_intx_mat*quads_mat[i])
            iota_mu_neg += Integer(periph_mat[1][i]*neg_intx_mat*quads_mat[i])

  
        pos_slopes.append((-iota_lambda_pos,iota_mu_pos))
        neg_slopes.append((-iota_lambda_neg,iota_mu_neg))
    return pos_slopes, neg_slopes

### this is almost the same as sign_bdy_matrix(), except we take the (entry-wise) absolute value
### of neg_intersection_mat(). This puts the orientation on quad edges described in Tillmann-- "Normal 
### surfaces in topologically finite 3-manifolds". The result are boundary slopes with signs indicating
### spinning direction
def spinning_bdy_maps(oriented_spun_surface, TN_wrapper):
    return signed_bdy_maps(oriented_spun_surface, TN_wrapper, True)

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

    M = TN_wrapper.manifold()
    peripheral_curve_matrices = TN_wrapper._peripheral_curve_mats

    intersection_mat = intersection_mat().apply_map(abs)
    slopes = [] 
    for c in range(T.countCusps()):
        periph_mat = peripheral_curve_matrices[c] # get the matrix that encodes meridian and longitude for cusp c.
        if quad_mat == None:
            quad_mat = oriented_quads_mat(s) # get the matrix that encodes quad coordinates for s.
        iota_lambda = QQ(sum([periph_mat[0][i]*intersection_mat*quad_mat[i] for i in range(T.size())]))
        iota_mu = QQ(sum([periph_mat[1][i]*intersection_mat*quad_mat[i] for i in range(T.size())]))
        slopes.append((-iota_lambda,iota_mu))
    return slopes


#def __boundary_slopes_(spun_surface, W=None):
#    s = spun_surface
#    T = s.triangulation()
#    if W == None:
#        M = snappy.Manifold(T.snapPea())
#        int_matrices = [intersection_mat(M, T, cusp) for cusp in range(T.countCusps())]
#    else:
#        M = W.manifold
#        int_matrices = W._peripheral_curve_mats
#    slopes = []
#    for cusp in range(T.countCusps()):
#        iota_lambda = 0
#        iota_mu = 0
#        mat = int_matrices[cusp]
#        for t in range(T.size()):
#            for v in range(4):
#                faces = [0,1,2,3]
#                faces.remove(v)
#                for f in faces:
#                    if 0 in [v,f]:
#                        j,k = sorted([v,f])
#                    else:
#                        verts = [0,1,2,3]
#                        verts.remove(v)
#                        verts.remove(f)
#                        j,k = verts
#                    iota_lambda += mat[0][t][4*v+f]*get_quads(s,t,j,k)
#                    iota_mu += mat[1][t][4*v+f]*get_quads(s,t,j,k)
#        slopes.append((-iota_lambda,iota_mu))
#    return slopes






