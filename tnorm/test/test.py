
import snappy
import tnorm.kernel.simplicial
import regina
import random


has_internal_homology = ['s789', 'v1539', 'v3209', 't06828', 't07788', 't08868', 't09692', 't12061','o9_09412', 'o9_23005', 'o9_23006', 'o9_31901', 'o9_32877', 'o9_36416', 'o9_36417', 'o9_36580', 'o9_40732', 'o9_41691', 'o9_41693', 'o9_44004', 'o9_44007', 'o9_44157']

closed_isosigs = ['kLLzvQQkacegjihihjjnkxmeqxjjcn', 'kLvPPPAkafedhhiigjjjqxnqxahnhj', 'kLvzQQzkccfehghgijjhnhqckcnahn', 'kLvzQQzkccfehghgijjhrhqckcnahr', 'kLLLAPPkcdgfehijijjhshaqpihkkk', 'kLLwvQQkacdgijhjiijjkalilmckwo', 'kLLzLPQkacegjihihjjjkxpuqxrwvv', 'kLLzLPQkacegjihihjjnkxpuqxrwvv', 'kLLzLPQkaceihjgijjijkxuetujjgs', 'kLLzLPQkaceihjgijjinkxuetujjgs']

links = ['L14n43208', 'L12n2026', 'L12n1054', 'L14n49678', 'L7a6', 'L8a18', 'L9a46']

  
nonfibered = ['K5a1','K6a3','K7a4','K7a5','K7a6','K7a3','K8a11','K8a18','K8a17','K8a10','K8a4','K8a9','K8a7','K8a1','K8a2','K9a27','K9a38','K9a35','K9a36','K9a23','K9a26','K9a8','K9a33','K9a39','K9a22','K9a34','K9a17','K9a10','K9a25','K9a24','K9a3','K9a21','K9a16','K9a4','K9a40','K9a18','K9a30','K9a32','K9a29','K9n5','K9n8','K10a75','K10a117','K10a113','K10a70','K10a65','K10a114','K10a64','K10a116','K10a43','K10a54','K10a33','K10a68','K10a115','K10a63','K10a108','K10a74','K10a60','K10a112','K10a57','K10a71','K10a61','K10a111','K10a58','K10a44','K10a34','K10a69','K10a55','K10a109','K10a19','K10a23','K10a5','K10a49','K10a29','K10a26','K10a30','K10a13','K10a82','K10a16','K10a80','K11a1','K11a2','K11a4','K11a6','K11a8','K11a10','K11a11','K11a12','K11a13','K11a16','K11a18','K11a20','K11a21','K11a23','K11a27','K11a29','K11a30','K11a31','K11a32','K11a36','K11a37','K11a38','K11a39','K11a41','K11a43','K11a45','K11a46','K11a48','K11a49','K11a50','K11a52','K11a54','K11a56','K11a58','K11a59','K11a60']


fibered = ['K3a1','K4a1','K5a2','K6a2','K6a1','K7a7','K7a2','K7a1','K8a8','K8a13','K8a6','K8a16','K8a3','K8a5','K8a15','K8a14','K8a12','K8n3','K8n1','K8n2','K9a41','K9a20','K9a14','K9a19','K9a2','K9a7','K9a15','K9a12','K9a5','K9a31','K9a1','K9a13','K9a6','K9a11','K9a28','K9a9','K9a37','K9n4','K9n3','K9n1','K9n2','K9n7','K9n6','K10a59','K10a56','K10a110','K10a107','K10a53','K10a35','K10a31','K10a52','K10a32','K10a25','K10a81','K10a15','K10a79','K10a2','K10a1','K10a41','K10a122','K10a38','K10a22','K10a10','K10a3','K10a27','K10a17','K10a78','K10a7','K10a83','K10a86','K10a11','K10a21','K10a106','K10a91','K10a24','K10a103','K10a104','K10a118','K10a72','K10a95','K10a66','K10a93','K10a100','K11a3','K11a5','K11a7','K11a9','K11a14','K11a15','K11a17','K11a19','K11a22','K11a24','K11a25','K11a26','K11a28','K11a33','K11a34','K11a35','K11a40','K11a42','K11a44','K11a47','K11a51','K11a53','K11a55','K11a57','K11a62','K11a66','K11a68','K11a71','K11a72','K11a73','K11a74','K11a76','K11a79','K11a80','K11a81','K11a82','K11a83','K11a86','K11a88','K11a92','K11a96','K11a99','K11a106','K11a108','K11a109','K11a112','K11a113']

fibered_knots = [random.choice(fibered) for i in range(10)]

nonfibered_knots = [random.choice(nonfibered) for i in range(10)]

tests = [closed_isosigs,has_internal_homology,links,fibered,nonfibered]

def test(do_tests=[0,1,2,3,4]):

    if 0 in do_tests:
        for name in closed_isosigs:
            print(name)
            T=regina.Triangulation3.fromIsoSig(name)
            W = tnorm.load(T, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
            for i in range(len(B.vertices())):
                assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'

            del3 = tnorm.kernel.simplicial.del3_matrix(W.triangulation())
            del2 = tnorm.kernel.simplicial.del2_matrix(W.triangulation())
            assert (del2*del3).is_zero(), W.manifold().name()+', force_simplicial=True'
            for i in range(W.qtons().size()):
                assert (del2*W.simplicial_class(i)).is_zero(), W.manifold().name()+', force_simplicial=True'
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=True'

    
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
            for i in range(len(B.vertices())):
                assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=False'

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
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=True'

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
            for i in range(len(B.vertices())):
                if B.vertices()[i].qtons_index() != None:
                    assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=False'

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
                if W.manifold().verify_hyperbolicity()[0]:
                    for i in range(W.qtons().size()):
                        assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=True'

    if 3 in do_tests:
        for name in fibered_knots:
            print(name)
            W = tnorm.load(name, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            assert not W.has_internal_homology(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
            assert B.polyhedron().is_compact(), W.manifold().name()+', force_simplicial=False'
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=False'

    if 4 in do_tests:
        for name in nonfibered_knots:
            print(name)
            W = tnorm.load(name, quiet=True)
            B = W.norm_ball
            assert B.polyhedron().dim() == W.triangulation().homologyH1().rank(), W.manifold().name()+', force_simplicial=False'
            assert not W.has_internal_homology(), W.manifold().name()+', force_simplicial=False'
            for i in range(W.qtons().size()):
                assert W.euler_char(i) == 2-2*W.genus(i)-W.num_boundary_comps(i), W.manifold().name()+', force_simplicial=False'
            assert B.polyhedron().is_compact(), W.manifold().name()+', force_simplicial=False'
            for i in range(len(B.vertices())):
                assert W.is_embedded(B.vertices()[i].qtons_index()) == True, W.manifold().name()+', force_simplicial=False'
            if W.manifold().verify_hyperbolicity()[0]:
                for i in range(W.qtons().size()):
                    assert W.euler_char(i) < 0, W.manifold().name()+', force_simplicial=False'        



