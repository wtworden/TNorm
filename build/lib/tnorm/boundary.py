
from tnorm.sage_types import *
from tnorm.regina_helpers import *
from tnorm.matrices import *

### count the number of boundary components of a spun normal surface
def num_boundary_comps_(spun_surface):
    numBoundaryComps = 0
    b = spun_surface.boundaryIntersections()
    for i in range(b.rows()):
        numBoundaryComps += gcd(b.entry(i,int(0)),b.entry(i,int(1)))
    return regina_to_sage_type(numBoundaryComps)


def get_quads(spun_surface,tet,v1,v2):
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_type(s.quads(tet,v2-int(1)))

def boundary_slopes_(spun_surface,M=None):
    s = spun_surface
    T = s.triangulation()
    if M == None:
        M = snappy.Manifold(T.snapPea())
    else:
        pass
    slopes = []
    for cusp in range(T.countCusps()):
        iota_lambda = 0
        iota_mu = 0
        mat = intersection_mat(M,T,cusp)
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







###################------------------------------------------------------------------


def get_surface_summary_homBdy(T,M,qns):
    angles = solve_lin_gluingEq(T)
    surface_types = []
    surface_types_w_index = []
    for i in range(qns.size()):
        st = Q_to_H1bdy(qns.surface(i)),euler_char_(qns.surface(i),angles)
        if st not in surface_types:
            surface_types.append(st)
            surface_types_w_index.append((i,st))
    return sorted(surface_types_w_index, key=lambda x: x[1][1].abs())

def get_surface_summary_normalBdy(T,M,qns):
    angles = solve_lin_gluingEq(T)
    surface_types = []
    surface_types_w_index = []
    for i in range(qns.size()):
        st = tuple([intersection_non_oriented(c,qns.surface(i),M) for c in range(T.countCusps())]),euler_char_(qns.surface(i),angles)
        if st not in surface_types:
            surface_types.append(st)
            surface_types_w_index.append((i,st))
    return sorted(surface_types_w_index, key=lambda x: x[1][1].abs())


