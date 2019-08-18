
import tnorm.constants
import tnorm.kernel.matrices
import tnorm.kernel.regina_helpers

intersection_mat = tnorm.kernel.matrices.intersection_mat
oriented_quads_mat = tnorm.kernel.matrices.oriented_quads_mat
regina_to_sage_type = tnorm.kernel.regina_helpers.regina_to_sage_type
ABS_SIGN_MATRIX = tnorm.constants.ABS_SIGN_MATRIX

from tnorm.sage_types import *



### count the number of boundary components of a spun normal surface
#def num_bdy_comps_unoriented_(spun_surface): 
#    numBoundaryComps = 0
#    b = spun_surface.boundaryIntersections()
#    for i in range(b.rows()):
#        numBoundaryComps += gcd(b.entry(i,int(0)),b.entry(i,int(1)))
#    return regina_to_sage_type(numBoundaryComps)


def get_quads(spun_surface,tet,v1,v2):
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_type(s.quads(tet,v2-int(1)))

### this map is almost the same as qtons_to_H1bdy(), except the SIGN_MATRIX is replaced by 
### SIGN_MATRIX.apply_map(abs). The result is the boundary slope for the embedded, possibly non-orientable,
### normal surface obtained by forgetting orientation. The reason this could be non-orientable is 
### that boundaries of the same slope that spin in oposite directions are cancelled. Thus the
### self-intersections coming from spinning are avoided, but at the cost of cutting and gluing
### in a way that does not respect orientation. This should return slopes which are (up to sign)
### the same as those returned by s.boundaryIntersections().
#def bdy_slopes_unoriented_(oriented_spun_surface, TN_wrapper=None):
#    s = oriented_spun_surface
#    T = s.triangulation()
#    if TN_wrapper == None:  # If we have already computed these thing, don't do it again.
#        M = snappy.Manifold(T.snapPea())
#        int_matrices = [intersection_mat(M, T, cusp) for cusp in range(T.countCusps())]
#    else:
#        M = TN_wrapper.manifold
#        int_matrices = TN_wrapper._intersection_matrices
#
#    sign_mat = ABS_SIGN_MATRIX # transverse orientation is ignored.
#    slopes = [] 
#    for c in range(T.countCusps()):
#        int_mat = int_matrices[c] # get the matrix that encodes meridian and longitude for cusp c.
#        quad_mat = oriented_quads_mat(s) # get the matrix that encodes quad coordinates for s.
#        iota_lambda = Rational(sum([int_mat[0][i]*sign_mat*quad_mat[i] for i in range(T.size())]))
#        iota_mu = Rational(sum([int_mat[1][i]*sign_mat*quad_mat[i] for i in range(T.size())]))
#        slopes.append((-iota_lambda,iota_mu))
#    return slopes


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






