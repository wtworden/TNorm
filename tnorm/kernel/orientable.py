from tnorm.utilities.regina_helpers import regina_to_sage_int
from tnorm.kernel.matrices import quads_mat

from tnorm.kernel.abstract_nbhd import abstract_nbhds

from tnorm.utilities.sage_types import Matrix


def orient(spun_normal_surface):
    S = spun_normal_surface
    T = S.triangulation()

    # We'll create a dictionary that keeps track of quad orientations wrt tets. For a given tet, 
    # with k quads, the quads are numbered from 0 to k-1, with the 0th quad being the one that is closest
    # to vertex 0 of tet. Then oriented_quads[tet] is a list of length k whose i^th entry is +1 if the 
    # orientation of the i^th quad agrees with the local orientation ("alignment"), and -1 if it disagrees.
    # At first, all entries are 0 to indicate that orientation is unknown (or really, not yet determined).


    abs_nbhds = abstract_nbhds(S)

    # abs_nbhds stores quads grouped according to the tetrahedron they are in. For below, it will be more
    # useful to store them according to which connected component of the intersection of the surface and 
    # the abstract neighborhood they are in.
    abstract_comps = []
    for e in abs_nbhds:
        N = abs_nbhds[e]
        for comp in N['components']:
            Qs = [i for i in range(len(comp)) if comp[i] != -1]
            abstract_comps.append({'tets':[N['tets'][i] for i in Qs], 'slopes':[N['slopes'][i] for i in Qs], 'align':[N['align'][i] for i in Qs], 'components':[comp[i] for i in Qs]})
    
    # iterator for assigning orientations to quads. See below.
    def iterate(comp, abstract_comps, oriented_quads, orientn_wrt_nbd=1):
        orientable = True

        # we iterate through the normal quadrilateral disks in the component comp, 
        # triangular disks are skips because we have already removed them from comp
        for i in range(len(comp['tets'])):

            # set orientation of this quad with respect to the tet
            orientn_wrt_tet = comp['align'][i]*orientn_wrt_nbd

            # if the orientation of the quad is not already set in "oriented_quads",
            # then set it to orientation found above
            if oriented_quads[comp['tets'][i]][comp['components'][i]] == 0:
                oriented_quads[comp['tets'][i]][comp['components'][i]] = orientn_wrt_tet
            
            # if there is already an orientation on the quad set in "oriented_quads",
            # then make sure it matches the orientation found above. If it doesn't,
            # then set orientalbe to False, and break.
            else:
                if oriented_quads[comp['tets'][i]][comp['components'][i]] != orientn_wrt_tet:
                    orientable = False
                    break

            # for each quad in comp, look for other components that it appears in, and add each 
            # of them, along with the orientation comp induces on them, to the list "new_comps"
            new_comps = []
            for other_comp in abstract_comps:
                for j in range(len(other_comp['tets'])):
                    if comp['tets'][i] == other_comp['tets'][j] and comp['components'][i] == other_comp['components'][j]:
                        new_orientn_wrt_nbd = orientn_wrt_tet*other_comp['align'][j]
                        new_comps.append((other_comp,new_orientn_wrt_nbd))
                        break

            # remove from abstract_comps all the components we just added to new_comps
            for new_comp in new_comps:
                abstract_comps.remove(new_comp[0])

            # for each component in new_comps, run iterate again. If new_comps is empty, then
            # we have exhausted this connected component, so we return
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
    oriented_quads_matrices = []
    while len(abstract_comps) > 0:
        # get a component. When we call iterate() below, we stay in this call to iterate
        # until we get orientations that are not compatible, or the component containing 
        # comp is completely oriented.
        oriented_quads = {}
        for i in range(T.size()):
            num_quads = regina_to_sage_int(S.quads(i,0) + S.quads(i,1) + S.quads(i,2))
            oriented_quads[i] = [0]*num_quads

        comp = abstract_comps.pop()
        abstract_comps, oriented_quads, orientable = iterate(comp, abstract_comps, oriented_quads)
        if not orientable:
            break

        # if orientable is still set to True, then this component is orientable, so we
        # compute the oriented_quads_mat for the component and store it in oriented_quads_matrices.
        if orientable:
            quads = quads_mat(S)
            emb_oriented_quads_mat = []
            for i in range(T.size()):
                row=[]
                for j in [0,1,2]:
                    if quads[i][j] != 0:
                        pos = sum([k for k in oriented_quads[i] if k==1])
                        neg = -sum([k for k in oriented_quads[i] if k==-1])
                        #assert pos+neg == quads[i][j]
                        row.append(pos)
                        row.append(neg)
                    else:
                        row.append(0); row.append(0)
                emb_oriented_quads_mat.append(row)
            oriented_quads_matrices.append(Matrix(emb_oriented_quads_mat))
    if orientable:
        return oriented_quads_matrices
    else:
        return False


def is_orientable(spun_normal_surface):
    S=spun_normal_surface
    if orient(S) != False:
        return True
    else:
        return False
