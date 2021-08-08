
from tnorm.utilities.sage_types import vector, Matrix
from tnorm.utilities.regina_helpers import regina_to_sage_int

"""
Let t be the number of tetrahedra for a triangulation T of M, and let v be the number of vertices (if M is cusped this is the
number of cusps). Then T has 4t/2=2t faces and which generate the dim-2 chain group C_2. A (possibly closed, not nec. embedded)
spun normal surface S is identified with an element of C_2, via an identification of quads with pairs of tetrahdera faces: an 
oriented quad q in a tetahedron t is isotopic to the pair of triangles of t (with outward tranverse orientation) that meet
at the edge that q misses, and which the transverse orientation of q points away from (see Segerman's "On spun normal and twisted
square surfaces"). By mapping all qtons to C2 via this scheme, we get a subspace K, which is the kernel of the boundary map 
from C2-->C1. We also get the image L of the boundary map C3-->C2 by taking the span of the images del3(t) for tetrahedra t in T.

Let dim(C2)=n, let dim(K)=k, and let dim(L)=l. let b_1,...,b_k be vectors in C2 forming a basis for K, and let c_1,...,c_l be 
vectors is C2 forming a basis for L. Let P be projection onto L. Then Q*b_i := b_i-P*b_i is in the orthogonal complement of L in K, as 
a subspace of C2. Let O=span{(Q*b_i)}_i. dim(O)=d should be the betti number of T, i.e., the dimension of H2 rel bdy.
Let {v_1,...,v_d} be a spanning set of vectors for O. Note that the v_i are not necesarily orthogonal, and there is no reason
to make them orthogonal, as we will map them to the stardand basis below.

Next find a map A:O-->QQ^d by solving A*B=Id_k where B=[v_1,...,v_d] is a matrix whose columns are the spanning vectors 
{v_1,...,v_d}. Now just map all normal surfaces down to H2 via A*Q. Note that the v_i come from cycles, so the basis elements 
for H2 are associated to normal surfaces.
"""


def standard_basis_vec(k, dim):
    return vector([0 if i != k else 1 for i in range(dim)])

def standard_basis_tup(k, dim):
    return tuple([0 if i != k else 1 for i in range(dim)])

def get_face_map_to_C2(T):
    tris = T.triangles()
    num_faces = len(tris)
    C2_face_map = {}
    for k in range(num_faces):
        t = tris[k]
        v = standard_basis_vec(k,num_faces)
        C2_face_map[(t.front().tetrahedron().index(),t.front().triangle())] = -v*(2*(t.front().triangle()%2)-1)
        C2_face_map[(t.back().tetrahedron().index(),t.back().triangle())] = v*(2*(t.front().triangle()%2)-1)

    return C2_face_map

def del3_matrix(T):
    rows = []
    C2_face_map = get_face_map_to_C2(T)
    for i in range(T.size()):
        row = vector([0 for _ in range(2*T.size())])
        for j in range(4):
            row += vector(C2_face_map[(i,j)])
        rows.append(row)
    mat = Matrix(rows).transpose()
    return mat


def get_quad_map_to_C2(T,C2_face_map):
    C2_quad_map = {}
    for i in range(T.size()):
        C2_quad_map[(i,0,False)] = C2_face_map[(i,0)]+C2_face_map[(i,1)]
        C2_quad_map[(i,0,True)] = C2_face_map[(i,2)]+C2_face_map[(i,3)]
        C2_quad_map[(i,1,False)] = C2_face_map[(i,0)]+C2_face_map[(i,2)]
        C2_quad_map[(i,1,True)] = C2_face_map[(i,1)]+C2_face_map[(i,3)]
        C2_quad_map[(i,2,False)] = C2_face_map[(i,0)]+C2_face_map[(i,3)]
        C2_quad_map[(i,2,True)] = C2_face_map[(i,1)]+C2_face_map[(i,2)]
    return C2_quad_map

#this maps a normal surface to C2
def map_qtons_to_C2(surface, C2_quad_map):
    T = surface.triangulation()
    image = vector([0 for i in range(2*T.size())])
    for i in range(T.size()):
        for j in range(3):
            for boo in [True,False]:
                image += regina_to_sage_int(surface.orientedQuads(i,j,boo))*C2_quad_map[(i,j,boo)]
    return image

# this returns vectors whose span is the subspace K of C2 in the notation above
def qtons_image_in_C2(TN_Wrapper, C2_quad_map): 
    W = TN_Wrapper
    image = []
    for i in range(W.qtons().size()):
        image.append(map_qtons_to_C2(W.qtons().surface(i),C2_quad_map))
    return image

def kernel_of_del2(TN_Wrapper, C2_quad_map):
    W = TN_Wrapper
    qtons_image = qtons_image_in_C2(W, C2_quad_map)
    return gram_schmidt(qtons_image)


# this returns vectors whose span is the subspace L of C2 in the notation above
def image_of_del3(TN_Wrapper, C2_face_map):
    W = TN_Wrapper
    T = W.triangulation()
    image = []
    for i in range(T.size()):
        tet = T.tetrahedron(i)
        bdy = 0
        for j in range(4):
            bdy += C2_face_map[(i,j)]
        image.append(bdy)
    return gram_schmidt(image)

# find qtons reps that realize longitudes.
def get_periph_basis_surfaces(TN_Wrapper):
    surfaces = []
    periph_basis_found = []
    W = TN_Wrapper
    if W.manifold().num_cusps() == 0:
        return []
    else:
        T = W.triangulation()
        num_cusps = W.manifold().num_cusps()
        betti = W.betti_number()
        periph_H2_std_vecs = [standard_basis_tup(k,num_cusps) for k in range(num_cusps)]
        for i in range(W.qtons().size()):
            H1bdy_slopes = W.H1bdy_slopes(i)
            periph_H2 = tuple([H1bdy_slopes[j][1] for j in range(len(H1bdy_slopes))])
            if periph_H2 in periph_H2_std_vecs:
                if periph_H2 not in periph_basis_found:
                    periph_basis_found.append(periph_H2)
                    surfaces.append(i)
        return surfaces

def get_nontrivial_closed(TN_Wrapper):
    W = TN_Wrapper
    return [i for i in range(W.qtons().size()) if W.qtons().surface(i).isCompact() and W.qtons().surface(i).cutAlong().isConnected()]


# returns vectors in C2 that form a basis for H2 (as a subspace of C2). This is O in the notation at the top.
def H2_as_subspace_of_C2(TN_Wrapper, C2_face_map, C2_quad_map):
    W = TN_Wrapper
    periph = get_periph_basis_surfaces(W)
    closed = get_nontrivial_closed(W)
    others = [i for i in range(W.qtons().size()) if ( not W.qtons().surface(i).isCompact() ) and ( i not in periph )]
    all_qtons = periph + closed + others
    qtons_image = qtons_image_in_C2(W,C2_quad_map)
    cycles = [qtons_image[i] for i in all_qtons]
    #K = gram_schmidt(cycles)
    L = image_of_del3(W, C2_face_map)

    proj = projection(L)
    non_orthog_H2_basis_in_C2 = find_spanning_subset([b-proj*b for b in cycles])

    return non_orthog_H2_basis_in_C2, proj, qtons_image

# same as gram_schmidt below, except we keep the original vectors instead of the orthogonalized vectors.
# So this just extracts the first d vectors in the list that span, where d=dim(span{vectors})
def find_spanning_subset(vectors): 
    dim = len(vectors[0])
    orthobasis = []
    spanning_subset = []
    for i in range(len(vectors)):
        orthobasis.append(vectors[i])
        spanning_subset.append(vectors[i])
        for j in range(0,len(orthobasis)-1):
            orthobasis[-1] -= ((vectors[i].dot_product(orthobasis[j]))/(orthobasis[j].dot_product(orthobasis[j])))*orthobasis[j]
        if orthobasis[-1].is_zero():
            orthobasis = orthobasis[:-1]
            spanning_subset = spanning_subset[:-1]
            
    return spanning_subset

# below implementation of gram-schmidt does not require vectors to be lin. ind. returns an orthogonal basis for 
# subspace generated by vectors. I.e., this returns span{vectors}.
def gram_schmidt(vectors): 
    dim = len(vectors[0])
    orthobasis = []
    for i in range(len(vectors)):
        orthobasis.append(vectors[i])
        for j in range(0,len(orthobasis)-1):
            orthobasis[-1] -= ((vectors[i].dot_product(orthobasis[j]))/(orthobasis[j].dot_product(orthobasis[j])))*orthobasis[j]
        if orthobasis[-1].is_zero():
            orthobasis = orthobasis[:-1]
            
    return orthobasis



def projection(vectors):
    A = Matrix(vectors).transpose()
    P = A*((A.transpose()*A).inverse())*A.transpose()
    return P


def del2_matrix(triangulation):
    bdy = []
    for tri in triangulation.triangles():
        row = [0] * triangulation.countEdges() # == countFaces(1)

        for i in range(3):
            edge_index = tri.edge(i).index()
            perm = tri.edgeMapping(i) # == tri.faceMapping(1,i)
            row[edge_index] += perm.sign()
        bdy.append(row)
    bdy = Matrix(bdy).transpose()
    return bdy




