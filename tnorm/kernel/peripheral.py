from tnorm.utilities.regina_helpers import regina_to_sage_int




### Pachner moves performed with SnapPy can result in a peripheral basis that has
### trivial loops, which throws off Euler characteristic calculations. So we need to
### check that each curve in the peripheral basis provided by SnapPy is connected. We
### also check that each meridian, longitude pair has algebraic intersection number +/- 1,
### which is a suffient condition to guarantee that they form a basis.

def alg_intersections(TN_wrapper):
    W = TN_wrapper
    P = W._peripheral_curve_mats
    Tri = W.triangulation()

    alg_int_numbers = []

    for cusp in range(W.num_cusps()):
        m_mat = P[cusp][0]
        l_mat = P[cusp][1]
        alg_int = 0

        for t in range(Tri.size()):
            row_t = m_mat[t]
            for v in range(4):
                for f in range(4):
                    if row_t[4*v+f] > 0:
                        b = row_t[4*v+f]

                        left, right = neighbors(v,f)
                        l_arcs = -(b + abs(row_t[4*v+left]) - abs(row_t[4*v+right]))/2
                        r_arcs = -(b + abs(row_t[4*v+right]) - abs(row_t[4*v+left]))/2

                        if l_mat[t][4*v+right] > 0:
                            alg_int += -(l_mat[t][4*v+right])*l_arcs/2


                        elif l_mat[t][4*v+right] < 0:
                            alg_int += -(l_mat[t][4*v+right])*l_arcs/2


                        if l_mat[t][4*v+left] > 0:
                            alg_int += (l_mat[t][4*v+left])*r_arcs/2


                        elif l_mat[t][4*v+left] < 0:
                            alg_int += (l_mat[t][4*v+left])*r_arcs/2

        alg_int_numbers.append(alg_int)
    
    return alg_int_numbers


def periph_basis_intersections(TN_wrapper):
    alg_ints = alg_intersections(TN_wrapper)
    good_basis = [abs(alg_ints[i]) == 1 for i in range(len(alg_ints))]
    if all(good_basis):
        return True, good_basis
    else:
        return False, 'peripheral basis for cusp {} does not have algebraic intersection +/- 1'.format(good_basis.index(False))


def periph_basis_connected(TN_wrapper):
    W = TN_wrapper
    P = W._peripheral_curve_mats
    Tri = W.triangulation()

    components = {}

    for cusp in range(W.num_cusps()):
        for curve in range(2):
            pmat = P[cusp][curve]
            
            visited = [(t,v) for t in range(Tri.size()) for v in range(4) if sum([abs(pmat[t][4*v+k]) for k in range(4)])>0]
            pos_weights = {(t,v,f):[1 for i in range(pmat[t][4*v+f])] for (t,v) in visited for f in range(4) if pmat[t][4*v+f] > 0}
            unvisited_arcs = sum([sum(value) for value in pos_weights.values()])

            t,v = visited[0]
            weights = [pmat[t][4*v+k] for k in range(4)]
            f = first_positive(weights)
            p = 1
            arc = (int(t),int(v),int(f),int(p))
            
            def iterate(component, arc, Tri):
                next_arc = arc_map(arc,pmat, Tri)
                component.append(next_arc)

                return component, next_arc


            component = [arc]
            unvisited_arcs -= 1
            pos_weights[(t,v,f)][p-1] -= 1
            component, arc = iterate(component, arc, Tri)

            while arc != component[0]:
                if unvisited_arcs == 0:
                    return False, 'peripheral curve {} on cusp {} not closed'.format(curve, cusp)
                unvisited_arcs -= 1
                pos_weights[(arc[0],arc[1],arc[2])][arc[3]-1] -= 1
                component, arc = iterate(component, arc, Tri)

            if unvisited_arcs != 0:
                return False, 'peripheral curve {} on cusp {} not connected'.format(curve, cusp)
            assert sum([sum([abs(w) for w in value]) for value in pos_weights.values()]) == 0
            components[cusp,curve] = component
    return True, components



def neighbors(v,f):
    faces = [1,2,3] if v==0 else [0,3,2] if v==1 else [0,1,3] if v==2 else [0,2,1]
    f_index = faces.index(f)
    left, right = faces[(f_index+1)%3], faces[f_index-1]
    return left, right

def first_positive(lst):
    for i in range(len(lst)):
        if lst[i] > 0:
            return i
    return None


def arc_map(arc, pmat, Tri):
    t, v, f, p = arc

    left, right = neighbors(v,f)
    etype = 1 if pmat[t][4*v + left] >= 0 else 2 if pmat[t][4*v + right] >= 0 else 3
    if etype == 1:
        ff = right
        p0 = pmat[t][4*v+left] + p
    elif etype == 2:
        ff = left
        p0 = p
    elif etype == 3:
        b = abs(pmat[t][4*v + left])
        c = abs(pmat[t][4*v + right])
        if p <= b:
            ff = left
            p0 = p
        elif p > b:
            ff = right
            p0 = p - b

    t0, v0, f0 = glued_to(t,v,ff, Tri)

    return (int(t0), int(v0), int(f0), int(p0))


def glued_to(t, v, f, Tri):
    tet = Tri.tetrahedron(t)
    t0 = tet.adjacentTetrahedron(f).index()
    v0 = tet.adjacentGluing(f)[v]
    f0 = tet.adjacentGluing(f)[f]

    return t0, v0, f0






