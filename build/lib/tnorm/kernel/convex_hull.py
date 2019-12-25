
import sympy as sp
from itertools import combinations

from sympy import ImmutableDenseMatrix as IMat
from sympy import Matrix as Mat
from sympy import Rational as Rational
from sympy import QQ


def get_test_points():
    points=[IMat([QQ(3),QQ(0),QQ(0)]),IMat([QQ(2),QQ(2),QQ(0)]),IMat([QQ(-2),QQ(2),QQ(0)]),IMat([QQ(3,2),QQ(-3,2),QQ(2)]),IMat([QQ(3,2),QQ(3,2),QQ(2)]),IMat([QQ(-3,2),QQ(-3,2),QQ(2)]),IMat([QQ(-3,2),QQ(3,2),QQ(2)]),IMat([QQ(0),QQ(0),QQ(3)])]
    all_points=[]
    for p in points:
        all_points.append(p)
        all_points.append(-p)
    return all_points

def get_test_face():
    p1 = Vertex(IMat([0,0,3]))
    p2 = Vertex(IMat([0,0,-3]))
    p3 = Vertex(IMat([3,0,0]))
    p4 = Vertex(IMat([-3,0,0]))
    e1 = Polyhedron(1,[p1,p3])
    e2 = Polyhedron(1,[p3,p2])
    e3 = Polyhedron(1,[p2,p4])
    e4 = Polyhedron(1,[p4,p1])
    f = Polyhedron(2,[e1,e2,e3,e4])
    return f

Finally, we have l A M = d I, so (l/d A) M = I and l/d A is the inverse
you seek.A N = d I

### finds the convex hull P of a set of points in d dimensions. Assumes the convex hull is a centrally symmetric polyhedron.
### 




def hull(points, dim):
    ## first we find a hyperplane B such that B cuts P in half along a maximal area cross-section.
    ## this allows us to just compute the convex hull of points to one side of B, then reflect.
    v_max = sorted([(ns(v),v) for v in points],key=lambda x:x[0])[-1][1]
    B_points = [v_max]
    #points.remove(v_max)
    #points.remove(-v_max)
    while len(B_points)<dim-1:
        A = Mat([v.transpose() for v in B_points]).transpose()
        Proj_P = A*(A.transpose()*A).inv()*A.transpose() #projection to the subspace generated by B_points
        perp_norms = dict([(v,ns(v-Proj_P*v)) for v in points])
        max_perp = sorted([(perp_norms[v],v) for v in points],key=lambda x:x[0])[-1][1]
        B_points.append(max_perp)
        #points.remove(max_perp)
        #points.remove(-max_perp)
    # remove all points to one side of B

    A = Mat([v.transpose() for v in B_points]).transpose()
    Proj_P = A*(A.transpose()*A).inv()*A.transpose()
    C = Mat([0 for i in range(dim-1)])
    A_t=A.transpose()

    solutions, param = A_t.gauss_jordan_solve(C)  # solve the system of equations, returns parametrization of solutions
    param_one = { tau:1 for tau in param } #set parameters to 0 to get a solution which will hopefully be an integer solution.
    normal = solutions.xreplace(param_one)
    to_remove = []
    for v in points:
        if normal.dot(v) < 0:
            to_remove.append(v)
    for v in to_remove:
        points.remove(v)
    # now we create our initial hull which includes points of B_points and their negatives, and the furthest point from the plane B.
    perp_norms = dict([(v,ns(v-Proj_P*v)) for v in points])
    max_perp = sorted([(perp_norms[v],v) for v in points],key=lambda x:x[0])[-1][1]
    B_points = [v for v in B_points]+[-v for v in B_points]
    F = hull(B_points,dim-1) # take the convex hull of points in B_points to get an initial face
    hull = double_cone_on_face(F,max_perp) # get the initial hull by taking the cone over F to max_perp
    outside = {}
    for f in hull.facets(dim-1):
        outside[f] = outside_of(f,points)

    return hull_points,points

class Polyhedron():
    def __init__(self,dim,facets=[]):
        self.dim = dim
        self._facets = {dim-1:facets}
        if self.dim > 1:
            for d in range(1,dim-1):
                self._facets[d]=[]
                for F in self._facets[dim-1]:
                    for f in F.facets(d):
                        if [v.coords() for v in f.vertices()] not in [[v.coords() for v in self._facets[d][i].vertices()] for i in range(len(self._facets[d]))]:
                            self._facets[d].append(f)
            verts=[]
            for F in self._facets[dim-1]:
                for v in F.vertices():
                    if v.coords() not in [w.coords() for w in verts]:
                        verts.append(v)
            self._facets[0] = tuple(verts)

    def __repr__(self):
        verts = ''
        for v in self.vertices():
            verts=verts+str(tuple(v.coords()[i] for i in range(len(v.coords()))))+','

        return "Convex hull of vertices: "+verts[:-1]

    def __hash__(self):
        return hash(self.vertices())

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and ( self.vertices() == other.vertices() )

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or ( self.vertices() != other.vertices() )

    def facets(self,d):
            return self._facets[d]

    def vertices(self):
            return self._facets[0]

    def is_subfacet(self,f1,f2):
        if f1 in f2.facets(f1.dim):
                return True
        else:
            return False

    def neighbors(self,vertex):
        neighbors = []
        for e in self.facets(1):
            if self.is_subfacet(vertex,e):
                for v in e.vertices():
                    if v != vertex:
                        neighbors.append(v)
        return neighbors


class Vertex():
    def __init__(self,coords):
        self._coords = coords
        self.dim = 0

    def __repr__(self):
        return str(tuple(self._coords[i] for i in range(len(self._coords))))

    def __hash__(self):
        return hash(self._coords)

    def __eq__(self,other):
        return ( self.__class__ == other.__class__ ) and self._coords == other._coords

    def __ne__(self,other):
        return ( self.__class__ != other.__class__ ) or self._coords != other._coords

    def coords(self,dim=0):
        return self._coords



def ns(vector):   ## norm squared
    return vector.norm()**2



def cone_on_face(face,vertex):
    if face.dim == 0:
        return Polyhedron(1,[face,vertex])
    else:
        facets = [face]
        for f in face.facets(face.dim-1):
            facets.append(cone_on_face(f,vertex))
        return Polyhedron(face.dim+1,facets)

def double_cone_on_face(face,vertex):
    if face.dim == 0:
        return Polyhedron(1,[face,vertex])
    else:
        facets = [face]
        for f in face.facets(face.dim-1):
            facets.append(cone_on_face(f,vertex))
        return Polyhedron(face.dim+1,facets)

def quick_hull(points, dim):
    # choose dim+1 points for the initial hull. We might as well pick the ones of largest norm.
    v_max = sorted([(ns(v),v) for v in points],key=lambda x:x[0])[-1][1]
    B_points = [v_max]
    while len(B_points)<dim:
        A = Mat([v.transpose() for v in B_points]).transpose()
        Proj_P = A*(A.transpose()*A).inv()*A.transpose() #projection to the subspace generated by B_points
        perp_norms = dict([(v,ns(v-Proj_P*v)) for v in points])
        max_perp = sorted([(perp_norms[v],v) for v in points],key=lambda x:x[0])[-1][1]
        B_points.append(max_perp)
    #O = Vertex(IMat([0 for i in range(dim)]))
    F = simplex(B_points,dim-1)
    facets = []
    for set1,set2 in partitions(F.vertices()):
        facets.append(simplex([-v.coords() for v in set1]+[v.coords() for v in set2],dim-1))
        facets.append(simplex([-v.coords() for v in set2]+[v.coords() for v in set1],dim-1))
    hull = Polyhedron(dim,facets)
    above = {} # this dictionary will associate with each face the vertices that are above it 
                 # (i.e., those having positive inner product with the outward normal of the face)
    for F in hull.facets(dim-1):
        above[F] = outside_of(F,points) 
    points = set([p for F in hull.facets(hull.dim-1) for p in above[F]]) # we can get rid of points not above at least one face
    while len(points) > 0: # as long as there are still points outside of the hull, expand faces to points that lie above them
        for F in hull.facets(dim-1):
            if len(above[F]) > 0:
                p = above[F][-1] # pick the point furthest from F (i.e., the last in the list)
        hull,points,above = expand(hull,p,points,above)

    #### ***** still need to combine co-planar faces. Compute a normal for each face, group by scalar multiple, combine faces in group.********

    return hull

def expand(hull,point,points,above):
    visible = [face for face in hull.facets(hull.dim-1) if point in above[face]] # faces visible from point
    horizons = [] # facets of dimension dim-2 that are on the boundary of the union of visible faces as seen from point
    for face in visible:
        for subface in face.facets(hull.dim-2):
            on_faces = [f for f in visible if hull.is_subfacet(subface,f)]
            if len(on_faces) == 1: # only add to horizons if it is not shared by multiple faces (so it's not in the interior)
                horizons.append(subface)
    faces = hull.facets(hull.dim-1) # all faces of current hull
    for face in visible:
        faces.remove(face) #remove visible faces as they are now covered by the new faces we'll create below
        del above[face] # also remove these from the above dictionary
    points.discard(point)
    for h in horizons:
        new_face = cone_on_face(h,Vertex(point)) # new faces come from coning over the horizons to point
        above[new_face] = outside_of(new_face,points) # populate above[new_face] with points that lie above new_face
        faces.append(new_face) # add to faces each new face

    hull = Polyhedron(hull.dim,faces) # construct the new hull
    points = set([p for F in hull.facets(hull.dim-1) for p in above[F]])
    return hull, points, above


def simplex(points,dim):
    assert len(points) == dim+1
    if dim == 1:
        P = Polyhedron(1,[Vertex(points[0]),Vertex(points[1])])
        return P
    else:
        facets = [simplex(f,dim-1) for f in list(combinations(points,dim))]
        P = Polyhedron(dim,facets)
        return P

def outside_of(face,points):
    outside = []
    basepoint = face.vertices()[0]  # pick a vertex to translate everything to the origin by
    pre_basis = face.neighbors(basepoint)
    dim = face.dim
    basis = [v.coords()-basepoint.coords() for v in pre_basis][:dim]
    A = Mat([v.transpose() for v in basis]).transpose()
    Proj_P = A*((A.transpose()*A).inv())*A.transpose()
    #C = Mat([0 for i in range(len(basis))])
    #A_t=A.transpose()

    #solutions, param = A_t.gauss_jordan_solve(C)  # solve the system of equations, returns parametrization of solutions
    #param_one = { tau:1 for tau in param } #set parameters to 0 to get a solution which will hopefully be an integer solution.
    #normal = solutions.xreplace(param_one)
    #if normal.dot(basepoint.coords()) < 0:
    #    normal = -normal
    normal = basepoint.coords() - Proj_P*basepoint.coords()
    for p in points:
        dot = normal.dot(p-basepoint.coords())
        if dot > 0:
            outside.append((dot,p))
    outside.sort(key=lambda x: x[0])
    return [pair[1] for pair in outside]  

class vector():
    pass

def partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        return [[ collection,[] ]]

    parts = []
    first = collection[0]
    for [set1, set2] in partitions(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        parts.append([set1,set2+[first]])
        parts.append([set1+[first],set2])
    return parts







