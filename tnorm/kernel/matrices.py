
from tnorm.sage_types import *
#from sympy import Matrix as Matrix
from tnorm.kernel.regina_helpers import get_oriented_quads, in_cusp


### for given cusp, return the matrix that gives intersection data for longitude and meridian.
### The Matrix mat_l has rows corresponding to tetrahedra, each row has length 16: the jth entry is 
### (j (mod 4) - 1) face of the the boundary triangle at the floor(j/4) vertex of the tet. A positive
### entry indicates that the longitude is going into the triangle, and a negative entry indicates
### it is exiting. Similarly for mat_m (meridian).
###                          v0    |     v1    |    v2     |    v3
###                    _                                               _
###               t0  | f0 f1 f2 f3 f0 f1 f2 f3 f0 f1 f2 f3 f0 f1 f2 f3 |
###               t1  |                       .                         |
###               .   |                       .                         |
###               .   |                       .                         |
###               .   |                       .                         |
###               .   |                       .                         |
###               tk  |_                                               _|

def peripheral_curve_mats(M,T,cusp):
    cusps = T.countCusps()
    tets = T.size()
    try:
        # for newer versions of SnapPy
        mat = M._get_cusp_indices_and_peripheral_curve_data()[1]
    except AttributeError:
        # for older versions of SnapPy
        mat = M._get_peripheral_curve_data()
    mat_l=Matrix([[mat[i][j]*int(in_cusp(T,i//4,j//4,cusp)) for j in range(16)] for i in range(2,4*tets,4)])
    mat_m=Matrix([[mat[i][j]*int(in_cusp(T,i//4,j//4,cusp)) for j in range(16)] for i in range(0,4*tets,4)])
    return (mat_l,mat_m)


### Quad coordinate of the surface as a Matrix. Row i corresponds to tetrahedron i, and each row has
### six entries corresponding to the six quad types: (0,1), (2,3), (0,2), (1,3), (0,3), (1,2) in that order.
def oriented_quads_mat(oriented_spun_surface):
    s = oriented_spun_surface
    T = s.triangulation()
    return Matrix([[get_oriented_quads(s,i,j,k) for (j,k) in [(0,1),(2,3),(0,2),(1,3),(0,3),(1,2)]] for i in range(T.size())])


### Same form as above matrix, but only picks up quads that have a boundary arc isotopic into the given
### cusp
def quad_bdy_mat(oriented_spun_surface, cusp):
    s = oriented_spun_surface
    T = s.triangulation()
    quads_mat = []
    for i in range(T.size()):
        row = []
        for (j,k) in [(0,1),(2,3),(0,2),(1,3),(0,3),(1,2)]:
            if ( in_cusp(T,i,j,cusp) or in_cusp(T,i,k,cusp) ):
                row.append(get_oriented_quads(s,i,j,k))
            else:
                row.append(0)
        quads_mat.append(row)
    return Matrix(quads_mat)


### Matrix that encodes how curves intersect with quad types, restricted to intersections
### that contribute to the positive boundary (see Cooper--Tillmann--Worden Fig 5).def pos_intx_matrix():
def pos_intersection_mat():
    return Matrix([[0, 0, 0, 0, 0, 0],
                   [1,0, 0, 0, 0, 0],
                   [0, 0, 1,0, 0, 0],
                   [0, 0, 0, 0, 1,0],
                   [1,0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0,0, 1],
                   [0, 0,0, 1, 0, 0],
                   [0, 0, 1,0, 0, 0],
                   [0, 0, 0, 0,0, 1],
                   [0, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1,0],
                   [0, 0,0, 1, 0, 0],
                   [0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0]])

### Matrix that encodes how curves intersect with quad types, restricted to intersections
### that contribute to the negative boundary (see Cooper--Tillmann--Worden Fig 5).
def neg_intersection_mat():
    return Matrix([[0, 0, 0, 0, 0, 0],
                   [0,-1, 0, 0, 0, 0],
                   [0, 0, 0,-1, 0, 0],
                   [0, 0, 0, 0, 0,-1],
                   [0,-1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0,-1, 0],
                   [0, 0,-1, 0, 0, 0],
                   [0, 0, 0,-1, 0, 0],
                   [0, 0, 0, 0,-1, 0],
                   [0, 0, 0, 0, 0, 0],
                   [-1, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0,-1],
                   [0, 0,-1, 0, 0, 0],
                   [-1, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0]])
### Matrix that encodes how curves intersect with quad types (see Cooper--Tillmann--Worden Fig 5).
def intersection_mat():
    return Matrix([[0, 0, 0, 0, 0, 0],
                   [1,-1, 0, 0, 0, 0],
                   [0, 0, 1,-1, 0, 0],
                   [0, 0, 0, 0, 1,-1],
                   [1,-1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0,-1, 1],
                   [0, 0,-1, 1, 0, 0],
                   [0, 0, 1,-1, 0, 0],
                   [0, 0, 0, 0,-1, 1],
                   [0, 0, 0, 0, 0, 0],
                   [-1, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1,-1],
                   [0, 0,-1, 1, 0, 0],
                   [-1, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0]])
def dot_product(v,w):
    return v.dot_product(w)

                  