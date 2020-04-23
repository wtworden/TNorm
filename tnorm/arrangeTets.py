
inf = float('inf')
from cmath import sqrt
import rhinoscriptsyntax as rs
import subprocess
import time
from os.path import expanduser
home = expanduser("~")
import pickle
import platform

#verts=[]
#for v in [v1,v2,v3]:
#    pt = rs.PointCoordinates(v)
#    pt = [pt[0],pt[1],pt[2]]
#    verts.append(pt)

    
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
        
    
def map_to_rs(point_on_ball):
    x = BH_inversion(point_on_ball)
    if x == inf:
        return inf
    else:
        return complex(x[0],x[1])
        
        
def map_to_ball(point_on_rs):
    x = point_on_rs
    if x != inf:
        x = [point_on_rs.real,point_on_rs.imag,0]
    print(x)
    return BH_inversion(x)
    
def map_to_standard(z,z1,z2,z3):
    w = ((z2-z3)*z-z1*(z2-z3))/((z2-z1)*z-z3(z2-z1))
    print(w)
    return w
    
def map_from_standard(z,z1,z2,z3):
    w = ((z3*(z1-z2))*z-(z1*(z3-z2)))/((z1-z2)*z+(z2-z3))
    print(w)
    return w

    
#Z1=map_to_rs(verts[0])
#Z2=map_to_rs(verts[1])
#Z3=map_to_rs(verts[2])
#print(Z1)

#W1 = complex(.5,sqrt(3)/2)
#W2 = complex(.5,-sqrt(3)/2)

#w1 = map_to_ball(map_from_standard(W1,Z1,Z2,Z3))
#w2 = map_to_ball(map_from_standard(W2,Z1,Z2,Z3))

#p1=rs.AddPoint(w1[0],w1[1],w1[2])
#p2=rs.AddPoint(w2[0],w2[1],w2[2])


#tet1 = rs.AddPoints([[w1[0],w1[1],w1[2]],v1,v2,v3])
#tet2 = rs.AddPoints([[w2[0],w2[1],w2[2]],v1,v2,v3])

if platform.system() == 'Darwin':
    filename = home + '/Library/Mobile Documents/com~apple~CloudDocs/3Ddrawings/temp/snappy_tmp.pickle'
    
elif platform.system() == 'Windows':
    filename = home + '\iCloudDrive\\3Ddrawings\\temp\\snappy_tmp.pickle'

with open(filename,'rb') as f:
    data = pickle.load(f)
    
shapes = data[0]
gluings = data[1]

