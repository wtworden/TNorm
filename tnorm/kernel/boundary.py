
from tnorm.kernel import matrices
from tnorm.utilities.regina_helpers import regina_to_sage_int

from tnorm.utilities.sage_types import Integer, QQ



def get_quads(spun_surface,tet,v1,v2):
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_int(s.quads(tet,v2-int(1)))


def signed_bdy_maps(oriented_spun_surface, peripheral_curve_mats, oriented_quads_mat=None, spinning=False):
    pass

### this maps an oriented spun normal surface to H_1(bdy M; \ZZ)
def inward_oriented_bdy(qtons, peripheral_curve_mats, oriented_quads_mat=None, spinning=False):
    s = qtons
    T = s.triangulation()
    quads_mat = oriented_quads_mat


    neg_intx_mat = matrices.NEG_INTERSECTION_MAT
    if spinning:
        neg_intx_mat = matrices.ABSNEG_INTERSECTION_MAT
    neg_slopes = []
    if quads_mat == None:
        quads_mat = matrices.oriented_quads_mat(s)

    neg_intx_times_quad_mat = (neg_intx_mat*quads_mat.transpose()).transpose()

    for c in range(T.countCusps()):
        periph_mat = peripheral_curve_mats[c] # get the matrix that encodes meridian and longitude for cusp c.
        
        iota_lambda_neg=0
        iota_mu_neg=0

        for i in range(T.size()):
            iota_lambda_neg += Integer(periph_mat[0][i]*neg_intx_times_quad_mat[i])
            iota_mu_neg += Integer(periph_mat[1][i]*neg_intx_times_quad_mat[i])

  
        neg_slopes.append((-iota_lambda_neg,iota_mu_neg))
    return neg_slopes


def outward_oriented_bdy(qtons, peripheral_curve_mats, oriented_quads_mat=None, spinning=False):
    s = qtons
    T = s.triangulation()
    quads_mat = oriented_quads_mat


    pos_intx_mat = matrices.POS_INTERSECTION_MAT
    pos_slopes = []                   # on a surface affects boundary orientations.
    if quads_mat == None:
        quads_mat = matrices.oriented_quads_mat(s)

    pos_intx_times_quad_mat = (pos_intx_mat*quads_mat.transpose()).transpose()

    for c in range(T.countCusps()):
        periph_mat = peripheral_curve_mats[c] # get the matrix that encodes meridian and longitude for cusp c.
        
        iota_lambda_pos=0
        iota_mu_pos=0

        for i in range(T.size()):
            iota_lambda_pos += Integer(periph_mat[0][i]*pos_intx_times_quad_mat[i])
            iota_mu_pos += Integer(periph_mat[1][i]*pos_intx_times_quad_mat[i])

  
        pos_slopes.append((-iota_lambda_pos,iota_mu_pos))
    return pos_slopes



### this is almost the same as sign_bdy_matrix(), except we take the (entry-wise) absolute value
### of neg_intersection_mat(). This puts the orientation on quad edges described in Tillmann-- "Normal 
### surfaces in topologically finite 3-manifolds". The result are boundary slopes with signs indicating
### spinning direction.
def outward_spinning_bdy(qtons, peripheral_curve_mats, oriented_quads_mat):
    return outward_oriented_bdy(qtons, peripheral_curve_mats, oriented_quads_mat, True)

def inward_spinning_bdy(qtons, peripheral_curve_mats, oriented_quads_mat):
    return inward_oriented_bdy(qtons, peripheral_curve_mats, oriented_quads_mat, True)


### this function computes spinning slopes for (unoriented) spun normal surfaces. 

def unoriented_spinning_slopes(spun_surface, peripheral_curve_mats, quads_mat=None):
    s = spun_surface
    T = s.triangulation()

    intx_mat = matrices.UNORIENTED_INTERSECTION_MAT
    slopes = []                   # on a surface affects boundary orientations.
    if quads_mat == None:
        quads_mat = matrices.quads_mat(s)

    intx_times_quad_mat = (intx_mat*quads_mat.transpose()).transpose()

    for c in range(T.countCusps()):
        periph_mat = peripheral_curve_mats[c] # get the matrix that encodes meridian and longitude for cusp c.
        
        iota_lambda=0
        iota_mu=0

        for i in range(T.size()):
            iota_lambda += Integer(periph_mat[0][i]*intx_times_quad_mat[i])
            iota_mu += Integer(periph_mat[1][i]*intx_times_quad_mat[i])

  
        slopes.append((-iota_lambda,iota_mu))
    return slopes












