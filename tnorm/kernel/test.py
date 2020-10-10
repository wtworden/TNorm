
import snappy
import tnorm.kernel.simplicial
import regina


has_internal_homology = ['s789', 'v1539', 'v3209', 't06828', 't07788', 't08868', 't09692', 't12061','o9_09412', 'o9_23005', 'o9_23006', 'o9_31901', 'o9_32877', 'o9_36416', 'o9_36417', 'o9_36580', 'o9_40732', 'o9_41691', 'o9_41693', 'o9_44004', 'o9_44007', 'o9_44157']

closed_isosigs = ['kLLzvQQkacegjihihjjnkxmeqxjjcn', 'kLvPPPAkafedhhiigjjjqxnqxahnhj', 'kLvzQQzkccfehghgijjhnhqckcnahn', 'kLvzQQzkccfehghgijjhrhqckcnahr', 'kLLLAPPkcdgfehijijjhshaqpihkkk', 'kLLwvQQkacdgijhjiijjkalilmckwo', 'kLLzLPQkacegjihihjjjkxpuqxrwvv', 'kLLzLPQkacegjihihjjnkxpuqxrwvv', 'kLLzLPQkaceihjgijjijkxuetujjgs', 'kLLzLPQkaceihjgijjinkxuetujjgs']

links = ['K10a110', 'L14n43208', 'L12n2026', 'K10n22', 'L12n1054', 'K10n13', 'L14n49678', 'L7a6', 'L8a18']

tests = [closed_isosigs,has_internal_homology,links]

def test(do_tests=[0,1,2]):

    if 0 in do_tests:
        for name in closed_isosigs:
            print(name)
            T=regina.Triangulation3.fromIsoSig(name)
            W = tnorm.load(T, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
#            assert B.polyhedron.is_compact(), W.manifold.name()+', force_simplicial=False'
        
            W = tnorm.load(name, quiet=True, force_simplicial_homology=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=True'
            del3 = tnorm.kernel.simplicial.del3_matrix(W.triangulation())
            del2 = tnorm.kernel.simplicial.del2_matrix(W.triangulation())
            assert (del2*del3).is_zero(), W.manifold().name()+', force_simplicial=True'
            for i in range(W.qtons().size()):
                assert (del2*W.simplicial_class(i)).is_zero(), W.manifold().name()+', force_simplicial=True'
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=True'
#            assert B.polyhedron.is_compact(), W.manifold.name()+', force_simplicial=True'
            for i in range(len(B.vertices())):
                assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'

    
    if 1 in do_tests:
        for name in has_internal_homology:
            print(name)
            W = tnorm.load(name, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            assert W.has_internal_homology(), W.manifold().name()+', force_simplicial=False'
            del3 = tnorm.kernel.simplicial.del3_matrix(W.triangulation())
            del2 = tnorm.kernel.simplicial.del2_matrix(W.triangulation())
            assert (del2*del3).is_zero(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert (del2*W.simplicial_class(i)).is_zero(), W.manifold().name()+', force_simplicial=False'
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
            assert B.polyhedron().is_compact(), W.manifold().name()+', force_simplicial=False'
        
            W = tnorm.load(name, quiet=True, force_simplicial_homology=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=True'
            assert W.has_internal_homology(), W.manifold().name()+', force_simplicial=True'
            del3 = tnorm.kernel.simplicial.del3_matrix(W.triangulation())
            del2 = tnorm.kernel.simplicial.del2_matrix(W.triangulation())
            assert (del2*del3).is_zero(), W.manifold().name()+', force_simplicial=True'
            for i in range(W.qtons().size()):
                assert (del2*W.simplicial_class(i)).is_zero(), W.manifold().name()+', force_simplicial=True'
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=True'
            assert B.polyhedron().is_compact(), W.manifold().name()+', force_simplicial=True'
            for i in range(len(B.vertices())):
                assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'
    
    if 2 in do_tests:
        for name in links:
            print(name)
            W = tnorm.load(name, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            assert not W.has_internal_homology(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
            assert B.polyhedron().is_compact(), W.manifold().name()+', force_simplicial=False'
        
            if W.manifold().num_cusps()>1:
                W_sim = tnorm.load(name, quiet=True, force_simplicial_homology=True)
                B_sim = W_sim.norm_ball
                assert B_sim.polyhedron().dim() == W_sim.triangulation().homologyH1().rank(), W_sim.manifold().name()+', force_simplicial=True'
                assert not W_sim.has_internal_homology(), W_sim.manifold().name()+', force_simplicial=True'
                del3 = tnorm.kernel.simplicial.del3_matrix(W_sim.triangulation())
                del2 = tnorm.kernel.simplicial.del2_matrix(W_sim.triangulation())
                assert (del2*del3).is_zero(), W_sim.manifold().name()+', force_simplicial=True'
                for i in range(W_sim.qtons().size()):
                    assert (del2*W_sim.simplicial_class(i)).is_zero(), W_sim.manifold().name()+', force_simplicial=True'
                    assert W_sim.euler_char(i) == 2-2*W_sim.genus(i)-W_sim.num_boundary_comps(i), W_sim.manifold().name()+', force_simplicial=True'
                assert B_sim.polyhedron().is_compact(), W_sim.manifold().name()+', force_simplicial=True'
                P = B.polyhedron()
                for i in range(len(B.vertices())):
                    assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'

                

