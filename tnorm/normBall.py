import TNorm

from tnorm.HashableDict import HashableDict
from tnorm.sage_types import *
from tnorm.hasse import get_hasse


class NormBall():
    def __init__(self,vertices,poly):
        self.verticesList = vertices
        self.polyhedron = poly
        self.numVertices = len(self.polyhedron.faces(0))
        self.numFaces = len(self.polyhedron.faces(self.polyhedron.dimension()-1))
        self.faceList = self.polyhedron.faces(self.polyhedron.dimension()-1)

    def vertex(self,key):    ### key can be index or coordinate of the vertex
        return self.verticesList.lookup(key)

    def show_hasse_diagram(self):
        G = get_hasse(self.polyhedron)
        G.show()

    def show_schlegel_diagram(self):
        G = Graphics()
        sp = self.polyhedron.schlegel_projection()
        p = sp.plot()
        G += p
        G.show(viewer='threejs', online=True)


    def hasse_diagram(self):
        return self.polyhedron.face_lattice().hasse_diagram()

    def plot(self):
        labels=Graphics()
        for i in range(self.numVertices):
            v=self.verticesList.lookup(i)
            labels+=text3d('S_{{{},{}}}'.format(v.genus,v.numPunctures),(v.coords[0],v.coords[1],v.coords[2]),color='black',dpi=200,fontsize='large',horizontal_alignment='center')
        p = self.polyhedron.plot(opacity=.7, edge_thickness=2)
        G = p+labels
        G.show(viewer='threejs',online=True)

class VerticesList():
    def __init__(self,vertices):

        self._verts = vertices
        self.by_coords = HashableDict(dict([(tuple(v.coords),v) for v in self._verts]))

    def lookup(self,key):
        if isinstance(key,tuple):
            return self.by_coords[key]
        elif isinstance(key,int) or isinstance(key, Integer):
            return self._verts[key]

    def detail(self):
        verts = []
        for v in self._verts:
            verts.append(v.__str__())
        return verts

    def __str__(self):
        string = ''
        for v in self._verts:
            string += v.__str__()+'\n' 
        return string


class Vertex():
    def __init__(self,index,surface_index,coords,num_punctures,euler,spinning_slopes,hom_slopes):
        self.index = index
        self.coords = coords
        self.surfaceIndex = surface_index
        self.numPunctures = num_punctures
        self.eulerChar = euler
        self.spinningSlopes = spinning_slopes
        self.H1BoundarySlopes = hom_slopes
        self.frac = (-1/self.eulerChar) if self.eulerChar != -1 else ''
        self.genus = (2- self.eulerChar - self.numPunctures)/2
    
    def __str__(self):
        return 'Vertex {}: represented by {}S_{},{} at {}, mapped from surface with index {}'.format(self.index,self.frac,self.genus,self.numPunctures,self.coords,self.surfaceIndex)

    def detail(self):
        return 'Vertex {}: represented by {}S_{},{} at {}, mapped from surface with index {}'.format(self.index,self.frac,self.genus,self.numPunctures,self.coords,self.surfaceIndex)

class Face():
    def __init__(self):
        pass