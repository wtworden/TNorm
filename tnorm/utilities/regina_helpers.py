
from .sage_types import Integer, QQ, Matrix



### return True if vertex "vert" of tetrahedron "tet" is in given cusp
def in_cusp(T,tet,vert,cusp):
    return T.tetrahedron(int(tet)).vertex(int(vert)) == T.cusp(int(cusp)).vertex()



### convert regina types (integers, rational, matrices) to sage types
def regina_to_sage_int(q):
    return Integer(q.longValue())

def regina_to_sage_rat(q):
    return QQ((q.numerator().longValue(),q.denominator().longValue()))

def regina_to_sage_mat(q):
    return Matrix([[regina_to_sage_int(q.entry(i,j)) for j in range(q.columns())] for i in range(q.rows())])



def get_oriented_quads(oriented_spun_surface,tet,v1,v2):
    assert v1 < v2
    s = oriented_spun_surface
    tet = int(tet)
    v2 = int(v2)
    if v1 == 0:
        return regina_to_sage_int(s.orientedQuads(tet,v2-int(1),True))
    elif v1 == 2:
        return regina_to_sage_int(s.orientedQuads(tet,int(0),False))
    elif v2 == 2:
        return regina_to_sage_int(s.orientedQuads(tet,int(2),False))
    else:
        return regina_to_sage_int(s.orientedQuads(tet,int(1),False))
    return "invalid input"

def get_quads(spun_surface,tet,v1,v2):
    assert v1 < v2
    assert v1 == 0
    s = spun_surface
    tet = int(tet)
    v2 = int(v2)
    return regina_to_sage_int(s.quads(tet,v2-int(1)))

