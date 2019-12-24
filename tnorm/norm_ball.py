from __future__ import print_function

import webbrowser
import time

from tnorm.x3d_to_html import make_x3d_html
from tnorm.sage_types import *
from tnorm.hasse import get_hasse
from tnorm.utilities import *


class NormBall():
    def __init__(self,vertices, polyhedron):
        self.vertices = vertices
        self.polyhedron = polyhedron
        self.dimension = polyhedron.dim()
        self.is_full_dimensional = self.polyhedron.is_full_dimensional()
        self.is_compact = self.polyhedron.is_compact()
        self.num_vertices = len(self.polyhedron.faces(0))
        self.num_faces = len(self.polyhedron.faces(self.polyhedron.dim()-1))
        self.face_list = self.polyhedron.faces(self.polyhedron.dim()-1)
        self._face_dict = {}
        self.confirmed = False

    ### sage does not index lines, rays, and vertices independently. We would like our indices to agree with those 
    ### of the sage Polyhedron, with lines and vertices indexed separately, so instead of using sages indices 
    ### via vertex.index(), we take the index of a vertex (or line) to be the index of the vertex (or line) in the tuple returned
    ### by self.polyhedron.vertices() (or self.polyhedron.lines()). This tuple is cached by sage, so the the order of the vertices in this 
    ### tuple should never change. The below function fetches this index.
    def index_of_poly_vert(self, v, is_line=False):   ## or line
        if is_line:
            verts = self.polyhedron.lines()
        else:
            verts = self.polyhedron.vertices()
        for i in range(len(verts)):
            if v == verts[i]:
                return i

    def vertex(self, index):
        return self.vertices[index]

    def show_schlegel_diagram(self):
        G = Graphics()
        sp = self.polyhedron.schlegel_projection()
        p = sp.plot()
        G += p
        G.show(viewer='threejs', online=True)

    @cached_property
    def get_hasse_diagram(self):
        G = get_hasse(self.polyhedron)
        return G

    def show_hasse_diagram(self):
        G = self.get_hasse_diagram
        G.show()

    def get_hasse_graph(self):
        return self.polyhedron.face_lattice().hasse_diagram()


    def _plot2d(self, dual=False):
        labels = Graphics()
        for i in range(self.num_vertices):
            v=self.vertices[i]
            if not dual:
                labels+=text('S_{},{}'.format(v.genus,v.num_boundary_comps), tuple(v.coords[j]*1.15 for j in range(len(v.coords))), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
            else:
                labels+=text('{}'.format(v.index), tuple(v.coords[j]*1.15 for j in range(len(v.coords))), color='black', dpi=200, fontsize='large', horizontal_alignment='center')
        plt = self.polyhedron.plot(dpi=200) + labels
        return plt

    def plot(self, viewer='x3d', online=False, opacity=1):
        dual = True if isinstance(self, DualNormBall) else False
        P = self.polyhedron
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
                for i in range(self.num_vertices):
                    v=self.vertices[i]
                    if not dual:
                        labels+=text3d('S_{},{}'.format(v.genus,v.num_boundary_comps),(v.coords[0]*1.15,v.coords[1]*1.15,v.coords[2]*1.15),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
                    else:
                        labels+=text3d('{}'.format(v.index),(v.coords[0]*1.15,v.coords[1]*1.15,v.coords[2]*1.15),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
                p = P.plot(opacity=opacity, edge_thickness=2)
                G = p+labels
                if with_dual:
                    pd = self.dual().plot(opacity=opacity, edge_thickness=2, color='blue')
                    G += pd
                G.show(viewer=viewer,online=online,axes=True)
        else:
            print('Can\'t plot polyhedron of dimension greater than 4.')

class TNormBall(NormBall):
    def __init__(self, vertices, rays, polyhedron):
        NormBall.__init__(self, vertices, polyhedron)
        self.rays = rays
        if not self.is_compact:
            self.polyhedron_mod_rays = Polyhedron(polyhedron.vertices_list())
        self.num_rays = len(self.polyhedron.rays())+2*len(self.polyhedron.lines())
        self.all_vertices_admissible = True
        for v in self.vertices:
             if not v.is_admissible:
                self.all_vertices_admissible = False
        
        for dim in range(1,self.polyhedron.dim()):
            self._face_dict[dim] = []
            for i in range(len(self.polyhedron.faces(dim))):
                face = self.polyhedron.faces(dim)[i]
                vertices = [self.vertex(self.index_of_poly_vert(v)) for v in face.vertices()]
                rays = [self.ray(2*self.index_of_poly_vert(l, True)+j) for l in face.lines() for j in [0,1]]
                self._face_dict[dim].append(Facet(i, dim, vertices, rays))



    def facets(self, dim):
        if dim == 0:
            return self.vertices
        else:
            return self._face_dict[dim]

    def facet(self, dim, index):
        if dim == 0:
            return self.vertices[index]
        else:
            return self._face_dict[dim][index]

    def num_facets(dim):
        return len(self.facets(dim))


    def ray(self, index):
        return self.rays[index]

    def project_to_norm(self):
        """ if the norm ball is not compact, project the infinite polyhedron to
            a subspace in which the Thurston norm is a true norm. While B.polyhedron_mod_rays
            is such a projection of P, it still lives in the same ambient space as P. This 
            function returns a polyhedron that is full dimensional in its ambient space (i.e.,
            the space is projected as well as the polyhedron). If B is compact, this just
            returns P (i.e., it does nothing). Also, if the dimension of P is < 4, this just returns
            B.polyhedron_mod_rays. This is becuase 
        """
        if B.polyhedron.dim() < 4:
            return self.polyhedron_mod_rays
        else:
            P = self.polyhedron_mod_rays
            subspace = P.change_ring(RR).ambient_space().subspace(P.vertices_list())
            basis = subspace.basis()
            projected_verts = [[vector(v).dot_product(b.normalized()) for b in basis] for v in P.vertices_list()]
            return Polyhedron(projected_verts)

class DualNormBall(NormBall):
    def __init__(self, vertices, polyhedron):
        NormBall.__init__(self, vertices, polyhedron)
        assert self.is_compact == True  # the dual norm ball should always be compact


class Vertex():
    def __init__(self, index, coords):
        self.index = index
        self.coords = coords


class AdornedVertex(Vertex):
    def __init__(self, index, surface_index, coords, num_boundary_comps, euler, boundary_slopes, is_ray, is_admissible, has_rep):
        Vertex.__init__(self, index, coords)
        self.surface_index = surface_index
        self.has_surface_rep = has_rep
        self.num_boundary_comps = num_boundary_comps 
        self.euler_char = euler
        self.boundary_slopes = boundary_slopes
        self._is_ray = is_ray
        if self.has_surface_rep:
            self.frac = (-1./self.euler_char) if self.euler_char != 0 else float('inf')
            self.genus = (2- self.euler_char - self.num_boundary_comps)/2
        else:
            self.frac = None
            self.genus = None
        self.is_admissible = is_admissible
    
    def __str__(self):
        frac = '' if self.frac == 1 else self.frac
        vert_or_ray = 'Ray' if self._is_ray else 'Vertex' 
        return '{} {}: represented by {}*S_{},{} at {}, mapped from surface with index {}'.format(vert_or_ray, self.index, self.frac, self.genus, self.num_boundary_comps, self.coords, self.surface_index)

    def __repr__(self):
        return self.__str__()

    def detail(self):
        return self.__str__()

class Ray(AdornedVertex):
    pass

class DualVertex(Vertex):
    def __init__(self, index, coords, dual_face):
        Vertex.__init__(self, index, coords)
        self.dual_face = dual_face


class Facet():
    def __init__(self, index, dim, vertices, rays):
        self.index = index
        self.dim = dim
        self.vertices = vertices
        self.rays = rays

    
    def __str__(self):
        return 'facet of dimension {} with vertices:\n{}'.format(self.dim,self.vertices)

    def __repr__(self):
        return self.__str__()






