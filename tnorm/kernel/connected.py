from tnorm.kernel.abstract_nbhd import abstract_nbhds





def is_connected(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    abs_nbhds = abstract_nbhds(S)

    comp_indices = [(k,j) for k in abs_nbhds for j in range(abs_nbhds[k]['components'].dimensions()[0])]

    stacks = {i:[] for i in range(T.size())}
    for j in abs_nbhds:
        comps = abs_nbhds[j]['components']
        for k in range(comps.dimensions()[1]):
            tet = abs_nbhds[j]['tets'][k]
            stacks[tet].append([(j,k),comps.column(k),abs_nbhds[j]['align'][k]])
            
    #return comp_indices, stacks