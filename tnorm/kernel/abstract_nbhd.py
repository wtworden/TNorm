from tnorm.utilities.regina_helpers import regina_to_sage_int


from tnorm.utilities.sage_types import Matrix



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

            # arranged quads is a vector of length width, whose components are in range(counts[ind]) for a quad, and -1 for a triangle.
            # e.g., vector [-1,0,1,-1] means in the tet at ind in cyclic ordering, top layer is a triangle, next down is
            # the 0th quad, next is 1st quad, next is a triangle (as we go along the edge from vertex 1 to vertex 0)
            arranged_quads = quad_range(counts[ind],align) + [-1]*(width-counts[ind])

            # where in the vertical ordering of the surface components are the first and last quad.
            first_quad, last_quad = 0,counts[ind]-1

            abstract_nbhd.append((tet,slope,align,arranged_quads,perm))

                
            # now continue, in the cyclic ordering, until we go all the way around the edge.
            while len(abstract_nbhd) < L:
                prev_slope = slope
                if counts[ind] != 0:
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
                    arranged_quads = [-1]*width
                else:
                    if prev_slope == 1 and slope == 1:
                        arranged_quads = pad_with_tris([-1]*(last_quad+1) + quad_range(counts[ind],align), width)
                    elif prev_slope == 1 and slope == -1:
                        arranged_quads = pad_with_tris([-1]*((prev_count + counts[ind]) + first_quad) + quad_range(counts[ind],align), width)
                    elif prev_slope == -1 and slope == 1:
                        arranged_quads = pad_with_tris([-1]*(first_quad) + quad_range(counts[ind],align), width)
                    elif prev_slope == -1 and slope == -1:
                        arranged_quads = pad_with_tris([-1]*(first_quad + counts[ind]) + quad_range(counts[ind],align), width)
                    first_quad, last_quad = first_non_neg(arranged_quads), last_non_neg(arranged_quads)
                abstract_nbhd.append((tet,slope,align,arranged_quads,perm))

        N = abstract_nbhd
        tets = [N[i][0] for i in range(len(N))]
        perms = [N[i][4] for i in range(len(N))]
        slopes =  [N[i][1] for i in range(len(N))]
        align = [N[i][2] for i in range(len(N))]
        components = Matrix([[N[i][3][j] for i in range(len(N))] for j in range(width)])

        # check to make sure in each component the quads alternate between positive and negative slopes. If
        # they don't then there is a problem, as in that case matching equations would not hold.
        for comp in components:
            total_sign=0
            for i in range(len(comp)):
                if comp[i] != -1:
                    total_sign += slopes[i]
                    assert total_sign in [-1,0,1], 'something is wrong---please let the developer know about this!'
            assert total_sign == 0, 'something is wrong---please let the developer know about this!'

            
        abst_nbhds[edge.index()] = {'perms':perms,'tets':tets,'slopes':slopes,'align':align,'components':components}

    return abst_nbhds








def quad_range(length,reversed):
    L = [i for i in range(abs(length))]
    if reversed == 1:
        L.reverse()
    return L


def first_non_neg(lst):
    i=0
    while lst[i] < 0 and i < len(lst)-1:
        i += 1
    if lst[i] >= 0:
        return i
    else:
        return None

def last_non_neg(lst):
    i = len(lst)-1
    while lst[i] < 0 and i > 0:
        i -= 1
    if lst[i] >= 0:
        return i
    else:
        return None



def pad_with_tris(lst, total_len):
    N = total_len - len(lst)
    return lst+[-1]*N

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

