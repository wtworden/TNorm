

import regina
from tnorm.sage_types import *


### identify which cusp the given vertex of the given tetrahedron is in
def get_cusp_index(T,tet,vertex):
    return T.tetrahedron(int(tet)).vertex(int(vertex)).boundaryComponent().index()

### return True if vertex "vert" of tetrahedron "tet" is in given cusp
def in_cusp(T,tet,vert,cusp):
    return get_cusp_index(T,tet,vert)==cusp


### convert regina types (integers, rational, matrices) to sage types
def regina_to_sage_type(q):
    if isinstance(q,regina.engine.Rational):
        return Rational(q.numerator().longValue(),q.denominator().longValue())
    elif isinstance(q,regina.engine.Integer):
        return Integer(q.longValue())
    elif isinstance(q,regina.engine.LargeInteger):
        return Integer(q.longValue())
    elif isinstance(q,regina.engine.MatrixInt):
        return Matrix([[regina_to_sage_type(q.entry(i,j)) for j in range(q.columns())] for i in range(q.rows())])



def get_oriented_quads(oriented_spun_surface,tet,v1,v2):
    assert v1 < v2
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
