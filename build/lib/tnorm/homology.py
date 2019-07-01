

from tnorm.sage_types import *
from tnorm.regina_helpers import *
from tnorm.matrices import *




def get_orientedQuads(oriented_spun_surface,tet,v1,v2):
    s = oriented_spun_surface
    tet = int(tet)
    v2 = int(v2)
    if v1 == 0:
        return regina_to_sage_type(s.orientedQuads(tet,v2-int(1),True))
    elif v1 == 2:
        return regina_to_sage_type(s.orientedQuads(tet,int(0),False))
    elif v2 == 2:
        return regina_to_sage_type(s.orientedQuads(tet,int(2),False))
    else:
        return regina_to_sage_type(s.orientedQuads(tet,int(1),False))
    return "invalid input"


### this is the old version of Q_to_H1bdy()
def intersection(cusp,oriented_spun_surface,M):
    s = oriented_spun_surface
    T = s.triangulation()
    T.size()
    mat = intersection_mat(M,T,cusp)
    print(mat)
    iota_lambda = 0
    iota_mu = 0
    for t in range(T.size()):
        for v in range(4):
            faces = [0,1,2,3]
            faces.remove(v)
            for f in faces:
                verts = [0,1,2,3]
                verts.remove(v)
                verts.remove(f)
                j,k = verts
                iota_lambda += mat[0][t][4*v+f]*get_orientedQuads(s,t,j,k)*(-1)
                iota_mu += mat[1][t][4*v+f]*get_orientedQuads(s,t,j,k)*(-1)
                j,k = sorted([v,f])
                iota_lambda += mat[0][t][4*v+f]*get_orientedQuads(s,t,j,k)
                iota_mu += mat[1][t][4*v+f]*get_orientedQuads(s,t,j,k)
    return (-iota_lambda,iota_mu)


### this maps an oriented spun normal surface to H_1(bdy M; \ZZ)
def Q_to_H1bdy(oriented_spun_surface,M=None):
    s = oriented_spun_surface
    T = s.triangulation()
    if M == None:
        M = snappy.Manifold(T.snapPea())
    else:
        pass
    s_mat = sign_matrix()
    slopes = []
    for c in range(T.countCusps()):
        int_mat = intersection_mat(M,T,c)
        quad_mat = orientedQuads_mat(s)
        iota_lambda = Rational(sum([int_mat[0][i]*s_mat*quad_mat[i] for i in range(T.size())]))
        iota_mu = Rational(sum([int_mat[1][i]*s_mat*quad_mat[i] for i in range(T.size())]))
        slopes.append((-iota_lambda,iota_mu))
    return slopes

### this maps an oriented spun normal surface to H_2(M,\partial M; \ZZ)
def homology_map(oriented_spun_surface,M=None):
    s = oriented_spun_surface
    if M == None:
        M = snappy.Manifold(s.triangulation().snapPea())
    else:
        pass
    slopes = Q_to_H1bdy(s,M)
    coords = vector([slopes[i][1] for i in range(len(slopes))])
    return coords


#################--------------------------------------------------------------------
