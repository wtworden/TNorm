

from tnorm.sage_types import *
from tnorm.regina_helpers import regina_to_sage_type


### compute the euler characteristic. If no angle structure is given, computes
### an angle structure to use for computation of Euler characteristic (i.e., the 
### second argument is optional). 
def euler_char_(spun_surface,angle_struct_matrix=None):
    T=spun_surface.triangulation()
    if angle_struct_matrix==None:
        angles = solve_lin_gluingEq(T)
    else:
        angles = angle_struct_matrix
    tets = T.size()
    ec = 0
    for i in range(tets):
        for j in range(3):
            ec -= angles[i][j]*regina_to_sage_type(spun_surface.quads(int(i),int(j)))
    return ec


##############-----------------------------------------------------------

### solve the linear gluing equations, including the linear part of the holonomy
### equations for meridians and longitudes. Returns a sage Matrix.
def solve_lin_gluingEq(T):
    GE = regina_to_sage_type(T.gluingEquations())
    for i in range(T.size()):
        # append equations log(z1)+log(z2)+log(z3)=\pi
        GE = GE.insert_row(GE.dimensions()[0],[Integer(0) for j in range(3*i)]+[Integer(1) for j in range(3)] + [Integer(0) for j in range(3*T.size()-3*(i+1))])
    edges = T.countEdges()
    cusps = T.countCusps()
    # 2pi for edge equations, 0 for cusp equations, pi for three dihedral angles of a tetrahedron (then divide out pi)
    v= vector([2 for i in range(edges)]+[0 for i in range(2*cusps)]+[1 for i in range(T.size())])
    assert v.length() == GE.dimensions()[0]
    angles_vec = GE.solve_right(v)  # solve the system of equations
    angles = Matrix([[angles_vec[j+i] for i in range(3)] for j in range(0,angles_vec.length(),3)])
    return angles











### It appears none of the below functions are being used.
##############-----------------------------------------------------------

def angle_struct_to_vec(regina_structure):
    tets = regina_structure.triangulation().size()
    return [regina.engine.Rational(regina_structure.angle(i,j)) for i in range(tets) for j in range(3)]

def eval_angle_struct(T,angle_struct):
    ge = T.gluingEquations()
    asv = angle_struct_to_vec(angle_struct)
    tets = T.size()
    edges = T.countEdges()
    cusps = T.countCusps()
    result = [regina.engine.Rational(0) for i in range(edges+2*cusps)]
    for i in range(edges+2*cusps):
        for j in range(3*tets):
            result[i] += regina.engine.Rational(ge.entry(i,j))*asv[j]
    return result

def find_good_structures(T,structures):
    edges = T.countEdges()
    cusps = T.countCusps()
    good_structures = []
    good = [regina.engine.Rational(2) for i in range(edges)]+[regina.engine.Rational(0) for i in range(2*cusps)]
    for i in range(structures.size()):
        if eval_angle_struct(T,structures.structure(i)) == good:
            good_structures.append(structures.structure(i))
    return good_structures

#############-------------------------------------------------------------------
