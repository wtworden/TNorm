from tnorm.kernel.abstract_nbhd import abstract_nbhds
from tnorm.kernel.matrices import quads_mat


from tnorm.utilities.sage_types import Graph

import regina


def connected_components(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    abst_nbhds = abstract_nbhds(S)
    quads_matrix = quads_mat(S)
    abst_nbhd_comps = [set([(abst_nbhds[k]['tets'][j],comp[j]) for j in range(len(comp)) if comp[j]!=-1]) for k in abst_nbhds for comp in abst_nbhds[k]['components']]

    anc_dict = dict([(i,abst_nbhd_comps[i]) for i in range(len(abst_nbhd_comps))])

    G = Graph(len(anc_dict))

    for i in G.vertices():
        for j in range(i+1, G.num_verts()):
            if len(anc_dict[i].intersection(anc_dict[j])) > 0:
                G.add_edge((i,j))

    components = []
    for C in G.connected_components():
        normal_coords = [int(0) for i in range(3*T.size())]
        quads = anc_dict[C[0]]
        for i in range(1,len(C)):
            quads = quads.union(anc_dict[C[i]])
        for (k,j) in quads:
            for i in range(3):
                if quads_matrix[k][i] != 0:
                    normal_coords[3*k+i] += 1
        F = regina.NormalSurface(T,regina.NS_QUAD, normal_coords)
        components.append(F)

    return components

def is_connected(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    abst_nbhds = abstract_nbhds(S)
    quads_matrix = quads_mat(S)
    abst_nbhd_comps = [set([(abst_nbhds[k]['tets'][j],comp[j]) for j in range(len(comp)) if comp[j]!=-1]) for k in abst_nbhds for comp in abst_nbhds[k]['components']]

    anc_dict = dict([(i,abst_nbhd_comps[i]) for i in range(len(abst_nbhd_comps))])

    G = Graph(len(anc_dict))

    for i in G.vertices():
        for j in range(i+1, G.num_verts()):
            if len(anc_dict[i].intersection(anc_dict[j])) > 0:
                G.add_edge((i,j))

    if G.is_connected():
        return True

    else:
        return False

def sum(R, S):

    T = R.triangulation()
    R_quads = quads_mat(R)
    S_quads = quads_mat(S)

    F_quads = R_quads + S_quads

    normal_coords = [int(F_quads[i][j]) for i in range(F_quads.dimensions()[0]) for j in range(3)]

    F = regina.NormalSurface(T, regina.NS_QUAD, normal_coords)

    return F

