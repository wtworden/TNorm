
from cmath import sqrt
import rhinoscriptsyntax as rs
import subprocess
import time
from os.path import expanduser
import pickle
import platform

home = expanduser("~")
inf = float('inf')


# Inversion that takes boundary of the ball model to boundary of upper half space model and vice versa   
def BH_inversion(point):
    x = point
    if x==[0,0,-1]:
        return inf
    elif x == inf:
        return [0,0,-1]
    else:
        ns = (x[0])**2+(x[1])**2+(x[2])**2
        print(ns)
        x[2] += 1
        print(x)
        ns = (x[0])**2+(x[1])**2+(x[2])**2
        print(ns)
        x=[2*c/ns for c in x]
        x[2] -= 1
        ns = (x[0])**2+(x[1])**2+(x[2])**2
        print(ns)
        return x
        

# map a point on S^2 to the Riemann sphere C U {infinity}
def map_to_rs(point_on_ball):
    x = BH_inversion(point_on_ball)
    if x == inf:
        return inf
    else:
        return complex(x[0],x[1])
        
# map a point of C U {\infty} = bdy(upper) to S^2 = bdy(ball)
def map_to_ball(point_on_rs):
    x = point_on_rs
    if x != inf:
        x = [point_on_rs.real,point_on_rs.imag,0]
    return BH_inversion(x)
    

# maps z1 to 1, z2 to 0, z3 to infinity    
def map_to_standard(z,z1,z2,z3):
    w = ((z1-z3)*z-z2*(z1-z3))/((z1-z2)*z-z3(z1-z2))
    return w

# maps 1 to z1, 0 to z2, infinity to z3        
def map_from_standard(z,z1,z2,z3):
    w = ((z3*(z2-z1))*z-(z2*(z3-z1)))/((z2-z1)*z+(z1-z3))
    return w



if platform.system() == 'Darwin':
    filename = home + '/Library/Mobile Documents/com~apple~CloudDocs/3Ddrawings/temp/snappy_tmp.pickle'
    
elif platform.system() == 'Windows':
    filename = home + '\\iCloudDrive\\3Ddrawings\\temp\\snappy_tmp.pickle'

with open(filename,'rb') as f:
    data = pickle.load(f)
    
shapes = [complex(shape[0],shape[1]) for shape in data[0]]
all_gluings = data[1]
num_tets = len(shapes)

verts=[]
for v in [v1,v2,v3]:
    pt = rs.PointCoordinates(v)
    pt = [pt[0],pt[1],pt[2]]
    verts.append(pt)
    
Z1=map_to_rs(verts[0])
Z2=map_to_rs(verts[1])
Z3=map_to_rs(verts[2])
#print(Z1)

W0 = shapes[init_tet]
#W2 = complex(.5,-sqrt(3)/2)

w0 = map_to_ball(map_from_standard(W0,Z1,Z2,Z3))
#w2 = map_to_ball(map_from_standard(W2,Z1,Z2,Z3))

p0=rs.AddPoint(w0[0],w0[1],w0[2])
#p2=rs.AddPoint(w2[0],w2[1],w2[2])

tet_verts = dict()
tet_verts[0] = [p0,v1,v2,v3]
tet0 = rs.AddPoints(tet_verts[0])
#tet2 = rs.AddPoints([[w2[0],w2[1],w2[2]],v1,v2,v3])

tets = []

for gluing in gluings:
    from_tet = int(gluing[0])
    to_tet = int(gluing[1])
    face = int(gluing[2])
    perm = dict([(i-4,int(gluing[i])) for i in range(4,8)])
    if face == 0:
        W = shapes[to_tet]
        i,j,k = 1,2,3
    elif face == 1:
        W = shapes[to_tet]
        W = (W-1)/W
        i,j,k = 0,2,3
    elif face == 2:
        W = shapes[to_tet]
        W = 1/(1-W)
        i,j,k = 0,1,3
    elif face == 3:
        W = shapes[to_tet]
        W = (W-1)/W
        i,j,k = 0,1,2
    v1,v2,v3 = tet_verts[from_tet][perm[i]], tet_verts[from_tet][perm[j]], tet_verts[from_tet][perm[k]]
    verts=[]
    for v in [v1,v2,v3]:
        pt = rs.PointCoordinates(v)
        pt = [pt[0],pt[1],pt[2]]
        verts.append(pt)
    Z1=map_to_rs(verts[0])
    Z2=map_to_rs(verts[1])
    Z3=map_to_rs(verts[2])
    w = map_to_ball(map_from_standard(W,Z1,Z2,Z3))
    p=rs.AddPoint(w[0],w[1],w[2])
    tet_verts[to_tet] = [p,v1,v2,v3]
    tet = rs.AddPoints(tet_verts[to_tet])
    tets.append(tet)



#class Tet:
#    def __init__(self, shape, gluings):
#        self.shape = shape
#        self.gluings = gluings
#
#    def glued_to(self,i):
#        return self.gluings[0][i], self.gluings[1][i]
#
#    def dihedral(i,j):
#        if (i,j) in [(0,1),(2,3)]:
#            return cmath.phase(self.shape)
#        elif (i,j) in [(1,2),(0,3)]:
#            return cmath.phase((self.shape-1)/(self.shape))
#        elif (i,j) in [(1,3),(0,2)]:
#            return cmath.phase(1/(1-self.shape))
#
#class Edge:
#    def __init__(self):
#        self.tets = []
#        self.angles = []
#
#class Face:
#    def __init__():
#        self.tets = []
#        self.is_glued = False
#
#class FunDom:
#    def __init__(self):
#        self.tets = []
#        self.edges = []
#        self.faces = []









