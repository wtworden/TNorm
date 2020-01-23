from __future__ import print_function
#from math import floor ## Warning: returns float

import sys
from tnorm.sage_types import *

def get_hasse(P):
    G = Graphics()
    FL = P.face_lattice()
    HD = FL.hasse_diagram()
    pos_dict, size_dict = get_vert_dicts(HD,P)
    HD.set_pos(pos_dict)
    half_dim = floor(P.dim()/2)
    w=len(P.faces(half_dim))
#    HD.set_latex_options(graphic_size=((P.dim()+3)*(w/4.),(P.dim()+3)*(w/6.)), vertex_size=1, edge_thickness=0.05, tkz_style="Shade",vertex_color="red")
    HD = HD.to_undirected()
    verts_by_dim = dict([(i,[v for v in HD.vertices() if v.dim()==i]) for i in range(-1,P.dim()+1)])
    min_dimension = min([i for i in verts_by_dim if len(verts_by_dim[i])!=0])
    for i in range(min_dimension,P.dim()+1):
        subgraph = HD.subgraph(verts_by_dim[i])
        subgraph.set_pos(dict([(v,pos_dict[v]) for v in verts_by_dim[i]]))
        size = size_dict[i]
        subgraph.relabel(relabel)
        if i==P.dim():
            p = subgraph.plot(dpi=200,vertex_size=200,vertex_labels=False,vertex_color="black")
        else:
            p = subgraph.plot(dpi=200,vertex_size=200*(size**1.7))
        G+=p
    HD.relabel(relabel)
    p = HD.plot(dpi=200,vertex_size=0,edge_thickness=1,vertex_labels=False)
    G+=p
    G.SHOW_OPTIONS['figsize'] = ((P.dim()+3)*(w/6.),(P.dim()+3)*(w/9.))
    G.axes(False)
    return G


def relabel(v):
    verts = v.vertices()
    length = len(verts)
    label = '{},'*(length-1)+'{}' if length>=2 else '{}'*length
#    label = '$\\text{{'+'{},'*(length-1)+'{}'+'}}$' if length>=2 else '$'+'{}'*length+'$'
    label = label.format(*[vert.index() for vert in verts])
    return label


def get_vert_dicts(G,P):
    dim = P.dim()
    verts_by_dim = dict([(i,[v for v in G.vertices() if v.dim()==i]) for i in range(-1,dim+1)])
    pos_dict = {}
    size_dict = {}
    levels = dim + 2
    levs = [1*i for i in range(levels)]
    v_top = verts_by_dim[dim][0]
    v_null = verts_by_dim[-1][0]
    pos_dict[v_top] = (.5,1-1./(dim+3))
    pos_dict[v_null] = (.5,1./(dim+3))
    size_dict[dim] = 1
    size_dict[-1] = 1
    half_dim = floor(dim/2.)
    g = 1./(len(verts_by_dim[half_dim]))
    for i in range(0,dim):
        verts = [v for v in G.vertices() if v in verts_by_dim[i]]
        verts.sort()
        maxlen = max([len(v.vertices()) for v in verts])
        pad = 2*abs(half_dim-i)*g + 2*g
        gap = (1-pad)/(len(verts)-1)
        if gap < g:
            gap = 1./len(verts)
            pad = gap
        size_dict[i] = maxlen
        for j in range(len(verts)):
            pos_dict[verts[j]] = (pad/2.+gap*j,(1./(dim+3)*(2+i)))
    return pos_dict, size_dict



#    [[1,1/3, 1/3, 1/3, 0],[1,1, 0, 0, 0],[1,0, 1/2, 0, 0],[1,0, 1/3, -1/3, -1/3],[1,0, 0, 1/2, 0],[1,0, 0, 0, 1],[1,0, 0, 0, -1],[1,0, 0, -1/2, 0],[1,0, -1/3, 1/3, 1/3],[1,0, -1/2, 0, 0],[1,-1, 0, 0, 0],[1,-1/3, -1/3, -1/3, 0]]

