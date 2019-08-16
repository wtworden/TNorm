from __future__ import print_function


import tnorm

from tnorm.sage_types import *
from tnorm.hasse import get_hasse
from tnorm.cached_property import *


class NormBall():
    def __init__(self,vertices,poly):
        self.vertices = vertices
        self.polyhedron = poly
        self.num_vertices = len(self.polyhedron.faces(0))
        self.num_faces = len(self.polyhedron.faces(self.polyhedron.dimension()-1))
        self.face_list = self.polyhedron.faces(self.polyhedron.dimension()-1)

    def vertex(self,index):
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

    def plot(self):
        labels=Graphics()
        for i in range(self.num_vertices):
            v=self.vertices[i]
            labels+=text3d('S_{{{},{}}}'.format(v.genus,v.num_punctures),(v.coords[0],v.coords[1],v.coords[2]),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
        p = self.polyhedron.plot(opacity=.7, edge_thickness=2)
        G = p+labels
        G.show(viewer='threejs',online=True)


class Vertex():
    def __init__(self,index,surface_index,coords,num_punctures,euler,qtons_slopes,hom_slopes):
        self.index = index
        self.coords = coords
        self.surface_index = surface_index
        self.num_punctures = num_punctures # this is for rep, not qtons. should have both.
        self.euler_char = euler
        self.qtons_slopes = qtons_slopes
        self.H1_boundary_slopes = hom_slopes
        self.frac = (-1/self.euler_char) if self.euler_char != -1 else ''
        self.genus = (2- self.euler_char - self.num_punctures)/2
    
    def __str__(self):
        return 'Vertex {}: represented by {}S_{},{} at {}, mapped from surface with index {}'.format(self.index,self.frac,self.genus,self.num_punctures,self.coords,self.surface_index)

    def detail(self):
        return 'Vertex {}: represented by {}S_{},{} at {}, mapped from surface with index {}'.format(self.index,self.frac,self.genus,self.num_punctures,self.coords,self.surface_index)

class Face():
    def __init__(self):
        pass




