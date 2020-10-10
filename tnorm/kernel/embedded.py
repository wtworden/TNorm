from tnorm.utilities import *
from tnorm.kernel.regina_helpers import *
from tnorm.kernel.matrices import oriented_quads_mat

def ends_embedded(qtons, TN_wrapper):
    W = TN_wrapper
    for i in range(W.num_cusps()):
        pos_bdy = W.spinning_slopes(qtons)[0][i]
        neg_bdy = W.spinning_slopes(qtons)[1][i]
        gcd_pos = gcd(pos_bdy[0],pos_bdy[1])
        gcd_neg = gcd(neg_bdy[0],neg_bdy[1])
        if gcd_pos != 0 and gcd_neg != 0:
            if not ( pos_bdy[0]/gcd_pos == neg_bdy[0]/gcd_neg and pos_bdy[1]/gcd_pos == neg_bdy[1]/gcd_neg ):
                return False
    return True

def is_admissible(qtons, TN_wrapper):
    W = TN_wrapper
    if W.allows_non_admissible():
        T = W.triangulation()
        for i in range(T.size()):
            if len([qtons.quads(i,j) for j in range(3) if qtons.quads(i,j) != 0]) > 1:
                return False
    else:
        return True

def abstract_nbhds(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    quad_counts = {}
    tets_around_edges = {}

    # First we'll gather, for each edge, information about the quads in the tetrahedra around the edge.
    # The tetrahedra are cylically ordered counter-clockwise w.r.t. the orientation of the edge with 
    # vertex on top. 
    for edge in T.faces(int(1)):

        # for each tet around edge, how many quads and what is slope (sign).
        quad_count = []

        # for each tet around edge, a tuple (tet_index, edge_type, edge embedding permutation)
        tets_around_edge = []

        embs = cyclically_order_embeddings(edge.embeddings(),T)
        for f in embs:
            t = f.tetrahedron()
            e = f.edge()
            p = f.vertices()
            tets_around_edge.append((t.index(),e,p))
            types = {0:[1,2], 5:[1,2], 1:[2,0], 4:[2,0], 2:[0,1], 3:[0,1]}
            for i in types[e]:
                count = regina_to_sage_int(S.quads(t.index(),i))
                if count != 0:
                    if i == types[e][1]:
                        count = -count
                    break
            quad_count.append(count)
        quad_counts[edge.index()] = quad_count
        tets_around_edges[edge.index()] = tets_around_edge

    # width is the number of connected components of the surface in the abstract nbhd at the edge,
    # ind is the index in the cyclic ordering around the edge where adding up crossing counts
    # from i to some later j gives the width. 
    widths = {}
    for edge in quad_counts:
        width, ind = 0, 0
        L = len(quad_counts[edge])
        for i in range(L):
            for j in range(i,L):
                this_width = int(sum(quad_counts[edge][i:j+1]))
                if abs(this_width) > width:
                    if this_width < 0:
                        width, ind = abs(this_width), j+1
                    else:
                        width, ind = this_width, i
        widths[edge] = (ind, width)

    # for each edge, abst_nbhds stores the info we need about the quads in its abstract nbhd.
    abst_nbhds = {}
    for edge in T.faces(int(1)):

        # for this edge, info we want to store about the abstract nbhd
        abstract_nbhd = []

        # get the counts for this edge, computed above
        counts = quad_counts[edge.index()]
        L = len(counts)

        # get ind and width for this edge, computed above
        ind, width = widths[edge.index()]

        # get tets_around_edge for this edge, computed above
        tets = tets_around_edges[edge.index()]

        # if there is anything in this abstract nbhd, continue
        if width > 0:

            # if there are no quads in the tet at ind in the ordering, then increasing ind until there are
            while counts[ind] == 0:
                ind = (ind + 1) % L
            assert counts[ind] > 0

            # For the first tet cyclically forward from ind containing quads, the quads will have positive
            # slope (because of how ind was chosen above), and will be as far toward vertex 1 of the edge as 
            # possible (i.e., in the top connected  component of the intersection of surface with abstract nbhd.)
            slope = 1
                
            tet, edge_type, perm = tets[ind]

            # For a permutation perm that described how edge is embedded in tet (or, equivalently, how tet 
            # is embedded in the abstract nbhd), this function tells us if the 0 vertex of tet is below the quad
            # (in which case alignment(perm)=-1) or above the quad (in which case alignment(perm)=1). This can
            # be interpreted as a transerve orientation of the quad, which gives a frame of reference. When we 
            # orient the quads, we can describe the orientation based on whether or not it agrees with this 
            # orientation.
            def alignment(perm):
                if perm[0]==0 or (slope==-1 and perm[3]==0) or (slope==1 and perm[2]==0):
                    return 1
                else:
                    return -1

            align = alignment(perm)

            # arranged quads is a vector of length width, whose components are 1 for a quad, and 0 for a triangle.
            # e.g., vector [0,1,1,0] means in the tet at ind in cyclic ordering, top layer is a triangle, next down is
            # a quad, next is a quad, next is a triangle (as we go along the edge from vertex 1 to vertex 0)
            arranged_quads = quad_range(counts[ind],align) + [None]*(width-counts[ind])

            # where in the vertical ordering of the surface components are the first and last quad.
            first_quad, last_quad = 0,counts[ind]-1

            abstract_nbhd.append((tet,slope,align,arranged_quads))

                
            # now continue, in the cyclic ordering, until we go all the way around the edge.
            while len(abstract_nbhd) < L:
                prev_slope = slope
                prev_count = counts[ind]
                ind = (ind + 1) % L

                if counts[ind] > 0:
                    slope = 1
                elif counts[ind] < 0:
                    slope = -1
                else: 
                    slope = prev_slope

                tet, edge_type, perm = tets[ind]
                align = alignment(perm)
                if counts[ind] == 0:
                    arranged_quads = [None]*width
                else:
                    if prev_slope == 1 and slope == 1:
                        arranged_quads = pad_with_None([None]*(last_quad+1) + quad_range(counts[ind],align), width)
                    elif prev_slope == 1 and slope == -1:
                        arranged_quads = pad_with_None([None]*((prev_count + counts[ind]) + first_quad) + quad_range(counts[ind],align), width)
                    elif prev_slope == -1 and slope == 1:
                        arranged_quads = pad_with_None([None]*(first_quad) + quad_range(counts[ind],align), width)
                    elif prev_slope == -1 and slope == -1:
                        arranged_quads = pad_with_None([None]*(first_quad + counts[ind]) + quad_range(counts[ind],align), width)
                    first_quad, last_quad = first_nonNone(arranged_quads), last_nonNone(arranged_quads)
                abstract_nbhd.append((tet,slope,align,arranged_quads))

        N = abstract_nbhd
        tets = [N[i][0] for i in range(len(N))]
        slopes =  [N[i][1] for i in range(len(N))]
        align = [N[i][2] for i in range(len(N))]
        components = [[N[i][3][j] for i in range(len(N))] for j in range(width)]
            
        abst_nbhds[edge.index()] = {'tets':tets,'slopes':slopes,'align':align,'components':components}

    return abst_nbhds

def quad_range(length,reversed):
    L = [i for i in range(abs(length))]
    if reversed == 1:
        L.reverse()
    return L


def orient(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    # We'll create a dictionary that keeps track of quad orientations wrt tets. For a given tet, 
    # with k quads, the quads are numbered from 0 to k-1, with the 0th quad being the one that is closest
    # to vertex 0 of tet. Then oriented_quads[tet] is a list of length k whose i^th entry is +1 if the 
    # orientation of the i^th quad agrees with the local orientation ("alignment"), and -1 if it disagrees.
    # At first, all entries are 0 to indicate that orientation is unknown (or really, not yet determined).
    oriented_quads = {}
    for i in range(T.size()):
        num_quads = regina_to_sage_int(S.quads(i,0) + S.quads(i,1) + S.quads(i,2))
        oriented_quads[i] = [0]*num_quads

    abs_nbhds = abstract_nbhds(S)

    # abs_nbhds stores quads grouped according to the tetrahedron they are in. For below, it will be more
    # useful to store them according to which connected component of the intersection of the surface and 
    # the abstract neighborhood they are in.
    abstract_comps = []
    for e in abs_nbhds:
        N = abs_nbhds[e]
        for j in range(len(N['components'])):
            Qs = [i for i in range(len(N['components'][j])) if N['components'][j][i] != None]
            abstract_comps.append({'tets':[N['tets'][i] for i in Qs], 'slopes':[N['slopes'][i] for i in Qs], 'align':[N['align'][i] for i in Qs], 'components':[N['components'][j][i] for i in Qs]})
    
    # iterator for assigning orientations to quads. See below.
    def iterate(comp, abstract_comps, oriented_quads, orientn_wrt_nbd=1):
        orientable = True
        for i in range(len(comp['tets'])):
            orientn_wrt_tet = comp['align'][i]*orientn_wrt_nbd
            if oriented_quads[comp['tets'][i]][comp['components'][i]] == 0:
                oriented_quads[comp['tets'][i]][comp['components'][i]] = orientn_wrt_tet
            else:
                if oriented_quads[comp['tets'][i]][comp['components'][i]] != orientn_wrt_tet:
                    orientable = False
                    break

            new_comps = []
            for other_comp in abstract_comps:
                for j in range(len(other_comp['tets'])):
                    if comp['tets'][i] == other_comp['tets'][j] and comp['components'][i] == other_comp['components'][j]:
                        new_orientn_wrt_nbd = orientn_wrt_tet*other_comp['align'][j]
                        new_comps.append((other_comp,new_orientn_wrt_nbd))
                        break

            for new_comp in new_comps:
                abstract_comps.remove(new_comp[0])

            for new_comp in new_comps:
                abstract_comps, oriented_quads, orientable = iterate(new_comp[0],abstract_comps,oriented_quads,new_comp[1])
                if not orientable:
                    return abstract_comps, oriented_quads, orientable
        return abstract_comps, oriented_quads, orientable

    # pick a first abstract component, and impose a transverse orientation. A transverse orientation is 
    # induced on each quad in the component. For each of these quads, we search for all remaining abstract
    # components in which they appear, and give that component the orientation induced by the orientation
    # on the quad. Iterate until we run out of components, or until we come to a component where the 
    # orientation induced by the quad that brought us there disagrees with the orientation of another quad
    # which has already been oriented (in which cases the surface is not orientable).
    while len(abstract_comps) > 0:
        comp = abstract_comps.pop()
        abstract_comps, oriented_quads, orientable = iterate(comp, abstract_comps, oriented_quads)
        if not orientable:
            break
    if orientable:
        return oriented_quads
    else:
        return False


def is_orientable(spun_normal_surface):
    if orient(S) != False:
        return True
    else:
        return False


def is_connected(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    abs_nbhds = abstract_nbhds(S)



def is_embedded(qtons,TN_wrapper):
    W = TN_wrapper
    S = qtons
    T = W.triangulation()

    if not W.manifold_is_closed():
        pos_bdy, neg_bdy = W.boundary_slopes(S)
        if not ends_embedded(S,W):
            return False

    oriented_projection = orient(S)

    if oriented_projection != False:
        q_mat = oriented_quads_mat(S)
        pos = [q_mat[i][0]+q_mat[i][2]+q_mat[i][4] for i in range(T.size())]
        pos_proj = [sum([j for j in oriented_projection[i] if j==1]) for i in range(T.size())]
        neg_proj = [sum([j for j in oriented_projection[i] if j==-1]) for i in range(T.size())]

        if pos==pos_proj or pos==neg_proj:
            return True
    return False


def first_nonNone(lst):
    i=0
    while lst[i] == None and i < len(lst)-1:
        i += 1
    if lst[i] != None:
        return i
    else:
        return None

def last_nonNone(lst):
    i = len(lst)-1
    while lst[i] == None and i > 0:
        i -= 1
    if lst[i] != None:
        return i
    else:
        return None



def pad_with_None(lst, total_len):
    N = total_len - len(lst)
    return lst+[None]*N

def cyclically_order_embeddings(edge_embeddings, T):
    """Order the edge embeddings cyclically
    """
    embs = edge_embeddings
    f = embs[0]
    e = f.edge()
    t = f.tetrahedron()
    dec = f.vertices()[0] < f.vertices()[1]
    ordered_embeddings = [f]
    emb_dict = dict([((f0.tetrahedron().index(),f0.edge(), f0.vertices()[0]<f0.vertices()[1]),f0) for f0 in embs])
    edge_dict = {(0,1):0, (0,2):1, (0,3):2, (1,2):3, (1,3):4, (2,3):5}
    gluing_dict = {(0,True):3, (0,False):2, (1,True):1, (1,False):3, (2,True):2, (2,False):1, (3,True):3, (3,False):0, (4,True):0, (4,False):2, (5,True):1, (5,False):0}
    while len(ordered_embeddings) < len(embs):
        gluing = t.adjacentGluing(gluing_dict[e,dec])
        t = t.adjacentTetrahedron(gluing_dict[e,dec])
        try:
            e = edge_dict[(gluing[f.vertices()[0]],gluing[f.vertices()[1]])]
            dec = True
        except KeyError:
            e = edge_dict[(gluing[f.vertices()[1]],gluing[f.vertices()[0]])]
            dec = False
        f = emb_dict[(t.index(),e,dec)]
        ordered_embeddings.append(f)
    return ordered_embeddings

















