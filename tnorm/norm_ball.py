from __future__ import print_function

import webbrowser
import time

from tnorm.x3d_to_html import make_x3d_html
from tnorm.sage_types import *
from tnorm.hasse import get_hasse
from tnorm.utilities import *


class NormBall(object):
    def __init__(self,vertices, polyhedron):
        self._vertices = vertices
        self._polyhedron = polyhedron
        self._dimension = polyhedron.dim()
        self._is_full_dimensional = self._polyhedron.is_full_dimensional()
        self._is_compact = self._polyhedron.is_compact()
        self._num_vertices = len(self._polyhedron.faces(0))
        self._num_faces = len(self._polyhedron.faces(self._polyhedron.dim()-1))
        self._face_list = self._polyhedron.faces(self._polyhedron.dim()-1)
        self._face_dict = {}
        self._confirmed = False

    def vertices(self):
        return self._vertices

    def polyhedron(self):
        return self._polyhedron

    def dim(self):
        return self._dimension

    def dimension(self):
        return self._dimension

    def is_full_dimensional(self):
        return self._is_full_dimensional

    def is_compact(self):
        return self._is_compact

    def num_vertices(self):
        return self._num_vertices

    def num_faces(self):
        return self._num_faces

    def face_list(self):
        return self._face_list

    def face_dict(self):
        return self._face_dict

    ### sage does not index lines, rays, and vertices independently. We would like our indices to agree with those 
    ### of the sage Polyhedron, with lines and vertices indexed separately, so instead of using sages indices 
    ### via vertex.index(), we take the index of a vertex (or line) to be the index of the vertex (or line) in the tuple returned
    ### by self.polyhedron.vertices() (or self.polyhedron.lines()). This tuple is cached by sage, so the the order of the vertices in this 
    ### tuple should never change. The below function fetches this index.
    def index_of_poly_vert(self, v, is_line=False):   ## or line
        if is_line:
            verts = self.polyhedron().lines()
        else:
            verts = self.polyhedron().vertices()
        for i in range(len(verts)):
            if v == verts[i]:
                return i

    def vertex(self, index):
        return self.vertices()[index]

    def show_schlegel_diagram(self):
        G = Graphics()
        sp = self.polyhedron().schlegel_projection()
        p = sp.plot()
        G += p
        G.show(viewer='threejs', online=True)

    @cached_property
    def get_hasse_diagram(self):
        G = get_hasse(self.polyhedron())
        return G

    def show_hasse_diagram(self):
        G = self.get_hasse_diagram
        G.show()

    def get_hasse_graph(self):
        return self.polyhedron().face_lattice().hasse_diagram()


    def _plot2d(self, dual=False):
        labels = Graphics()
        for i in range(self.num_vertices()):
            v=self.vertices()[i]
            if not dual:
                if self.dimension() ==1:
                    labels+=text('S_{},{}'.format(v.genus(),v.num_boundary_comps()), (v.coords()[0]*1.15,0), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
                else:
                    labels+=text('S_{},{}'.format(v.genus(),v.num_boundary_comps()), tuple(v.coords()[j]*1.15 for j in range(len(v.coords()))), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
            else:
                if self.dimension() ==1:
                    labels+=text('{}'.format(v.index()), (v.coords()[0]*1.15,0), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
                else:
                    labels+=text('{}'.format(v.index()), tuple(v.coords()[j]*1.15 for j in range(len(v.coords()))), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
        plt = self.polyhedron().plot(dpi=200) + labels
        return plt

    def plot(self, viewer='x3d', online=False, opacity=1):
        dual = True if isinstance(self, DualNormBall) else False
        P = self.polyhedron()
        if not P.is_compact():
            P = Polyhedron(P.vertices_list()) # we can't plot an infinite polyhedron, so just plot the finite part.
        if P.dim() in [1,2]:
            plt = self._plot2d(dual)
            plt.show()
        elif P.dim() in [3,4]:
            if viewer == 'x3d':
                with make_temp_directory() as temp_dir:
                    filename = make_x3d_html(self, online=online, directory=temp_dir, dual=dual)
                    webbrowser.open('file://'+filename)
                    time.sleep(5)
            else:
                labels=Graphics()
                for i in range(self.num_vertices()):
                    v=self.vertices()[i]
                    if not dual:
                        labels+=text3d('S_{},{}'.format(v.genus(),v.num_boundary_comps()),(v.coords()[0]*1.15,v.coords()[1]*1.15,v.coords()[2]*1.15),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
                    else:
                        labels+=text3d('{}'.format(v.index()),(v.coords()[0]*1.15,v.coords()[1]*1.15,v.coords()[2]*1.15),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
                p = P.plot(opacity=opacity, edge_thickness=2)
                G = p+labels
                G.show(viewer=viewer,online=online,axes=True)
        else:
            print('Can\'t plot polyhedron of dimension greater than 4.')

class TNormBall(NormBall):
    def __init__(self, vertices, rays, polyhedron):
        NormBall.__init__(self, vertices, polyhedron)
        self._rays = rays
        if not self.is_compact():
            self._polyhedron_mod_rays = Polyhedron(polyhedron.vertices_list())
        else:
            self._polyhedron_mod_rays = polyhedron
        self._num_rays = len(self.polyhedron().rays())+2*len(self.polyhedron().lines())
        self._all_vertices_admissible = True
        for v in self.vertices():
             if not v.is_admissible():
                self._all_vertices_admissible = False
        
        for dim in range(1,self.polyhedron().dim()):
            self._face_dict[dim] = []
            for i in range(len(self.polyhedron().faces(dim))):
                face = self.polyhedron().faces(dim)[i]
                vertices = [self.vertex(self.index_of_poly_vert(v)) for v in face.vertices()]
                rays = [self.ray(2*self.index_of_poly_vert(l, True)+j) for l in face.lines() for j in [0,1]]
                self._face_dict[dim].append(Facet(i, dim, vertices, rays))

    def polyhedron_mod_rays(self):
        return self._polyhedron_mod_rays

    def num_rays(self):
        return self._num_rays

    def all_vertices_admissible(self):
        return self._all_vertices_admissible

    def simplicial_class_of_vertex(self,vertex_index):
        return self.vertices[vertex_index].simplicial_class

    def facets(self, dim):
        if dim == 0:
            return self.vertices()
        else:
            return self._face_dict[dim]

    def facet(self, dim, index):
        if dim == 0:
            return self.vertices()[index]
        else:
            return self._face_dict[dim][index]

    def num_facets(dim):
        return len(self.facets(dim))

    def rays(self):
        return self._rays

    def ray(self, index):
        return self._rays[index]

    def project_to_norm(self):
        """ if the norm ball is not compact, project the infinite polyhedron to
            a subspace in which the Thurston norm is a true norm. While B.polyhedron_mod_rays
            is such a projection of P, it still lives in the same ambient space as P. This 
            function returns a polyhedron that is full dimensional in its ambient space (i.e.,
            the space is projected as well as the polyhedron). If B is compact, this just
            returns P (i.e., it does nothing). Also, if the dimension of P is < 4, this just returns
            B.polyhedron_mod_rays. This is becuase 
        """
        if self.polyhedron().dim() < 4:
            return self.polyhedron_mod_rays()
        else:
            P = self.polyhedron_mod_rays()
            subspace = P.change_ring(RR).ambient_space().subspace(P.vertices_list())
            basis = subspace.basis()
            projected_verts = [[vector(v).dot_product(b.normalized()) for b in basis] for v in P.vertices_list()]
            return Polyhedron(projected_verts)

class DualNormBall(NormBall):
    def __init__(self, vertices, polyhedron):
        NormBall.__init__(self, vertices, polyhedron)
        assert self.is_compact() == True  # the dual norm ball should always be compact


class Vertex(object):
    def __init__(self, index, coords):
        self._index = index
        self._coords = coords

    def coords(self):
        return self._coords

    def index(self):
        return self._index


class NBVertex(Vertex):
    def __init__(self, index, qtons_index, coords, TN_wrapper, norm_minimizer_info=None):
        Vertex.__init__(self, index, coords)

        # the qtons surface that is a multiple of this vertex.
        self._qtons_index = qtons_index

        self._TN_wrapper = TN_wrapper

        if qtons_index != None:
            self._has_qtons_rep = True
            self._norm_minimizer_bdy_comps = self._TN_wrapper.num_boundary_comps(self._qtons_index)
            self._norm_minimizer_euler = self._TN_wrapper.euler_char(self._qtons_index)
            self._norm_minimizer_bdy_slopes = self._TN_wrapper.boundary_slopes(self._qtons_index)
            self._norm_minimizer_genus = self._TN_wrapper.genus(self._qtons_index)  
            self._norm_minimizer_is_embedded = self._TN_wrapper.is_embedded(self._qtons_index)
            self._simplicial_class = self._TN_wrapper.simplicial_class(self._qtons_index)
            self._is_admissible = self._TN_wrapper.is_admissible(self._qtons_index)    
        else:
            self._has_qtons_rep = False

            # norm_minimizer_info is a tuple containing (num_boundary_comps, euler_char, bdy_slopes).
            # Since there is no spun normal surface representing the norm minimizer, this info does not
            # correspond to a surface in TN_wrapper.qtons(), but does correspond to the norm minimizer
            # for the primitive integer point lying over the vertex.
            self._norm_minimizer_bdy_comps = norm_minimizer_info[0]
            self._norm_minimizer_euler = norm_minimizer_info[1]
            self._norm_minimizer_bdy_slopes = norm_minimizer_info[2]
            self._norm_minimizer_genus = (2-self._norm_minimizer_euler-self._norm_minimizer_bdy_comps)/2
            self._simplicial_class = None
            self._is_admissible = None
            self._norm_minimizer_is_embedded = None

        self._is_ray = False
        self._factor = QQ(-1/self.euler_char()) if self.euler_char() != 0 else float('inf')

        # the index of the spun normal surface that is norm minimizing and is the primitive integer point on the ray through this vertex.
    def qtons_index(self):
        return self._qtons_index

    def simplicial_class(self):
        """return the simplicial class of the norm minimizing surface that is a multiple of this vertex 
            (i.e., it is a primitive integer point lying on the ray over this vertex)
        """
        return self._simplicial_class

    def TN_wrapper(self):
        return self._TN_wrapper

    def is_ray(self):
        return self._is_ray

    def genus(self):
        """return the genus of the norm minimizing surface that is a multiple of this vertex (i.e., it 
            is a primitive integer point lying on the ray over this vertex)
        """
        return self._norm_minimizer_genus

    def num_boundary_comps(self):
        """return the number of boundary components of the norm minimizing surface that is a multiple of 
            this vertex (i.e., it is a primitive integer point lying on the ray over this vertex)
        """
        return self._norm_minimizer_bdy_comps

    def euler_char(self):
        """return the Euler characteristic of the norm minimizing surface that is a multiple of this vertex 
            (i.e., it is a primitive integer point lying on the ray over this vertex)
        """
        return self._norm_minimizer_euler

    def boundary_slopes(self):
        """return the boundary slopes of the norm minimizing surface that is a multiple of this vertex 
            (i.e., it is a primitive integer point lying on the ray over this vertex)
        """
        return self._norm_minimizer_bdy_slopes   

    def is_embedded(self):
        return self._norm_minimizer_is_embedded

    def is_admissible(self):
        return self._is_admissible
    
    def __str__(self):
        frac = '' if self._factor == 1 else self._factor
        vert_or_ray = 'Ray' if self._is_ray else 'Vertex' 
        return '{} {}: represented by ({})*S_{},{} at {}, mapped from surface with index {}'.format(vert_or_ray, self.index(), frac, self.genus(), self.num_boundary_comps(), self.coords(), self.qtons_index())

    def __repr__(self):
        return self.__str__()

    def detail(self):
        return self.__str__()

class Ray(NBVertex):
    def __init__(self, index, qtons_index, coords, TN_wrapper):
        NBVertex.__init__(self, index, qtons_index, coords, TN_wrapper)
        self._is_ray = True

class DualVertex(Vertex):
    def __init__(self, index, coords, dual_face):
        Vertex.__init__(self, index, coords)
        self._dual_face = dual_face

    def dual_face(self):
        return self._dual_face

class Facet(object):
    def __init__(self, index, dim, vertices, rays):
        self._index = index
        self._dim = dim
        self._vertices = vertices
        self._rays = rays

    def index(self):
        return self._index

    def dim(self):
        return self._dim

    def vertices(self):
        return self._vertices

    def rays(self):
        return self._rays

    
    def __str__(self):
        return 'facet of dimension {} with vertices:\n{}'.format(self.dim(),self.vertices())

    def __repr__(self):
        return self.__str__()






