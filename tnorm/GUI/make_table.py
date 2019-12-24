from __future__ import print_function
import tnorm

import tkinter as tk

from tkinter import ttk
from math import ceil



def make_qtons_table(tnorm_app):
#    try:
    qd = tnorm_app.wrapper.qtons_info
    treedata = [(i, qd[i]['image_in_H2'], qd[i]['euler_char'], str(qd[i]['boundary_slopes']), qd[i]['num_boundary_comps'], qd[i]['genus'], qd[i]['is_norm_minimizing'], qd[i]['over_face'], str(qd[i]['spinning_slopes'])) for i in qd]
    column_names = ('qtons index', 'image in H2', 'euler char', 'boundary slopes', 'num bdy comps', 'genus','norm minimizing?','over face','spinning_slopes')
    make_treeview_table(tnorm_app.QtonsTab,treedata,column_names, True)
#    except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()


def make_all_vertices_table(tnorm_app, frame):
    #try:
    make_vertices_table(tnorm_app, frame, 'all', True)
    #except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()

def make_vertices_table(tnorm_app, frame, vert_list, with_scrollbar=True):
    W = tnorm_app.wrapper
    B = tnorm_app.ball
    if vert_list == 'all':
        vert_list = range(len(B.vertices))
    height = len(vert_list)
    width = 9
    verts = {i:B.vertices[i] for i in vert_list}
    column_names = ('vertex','qtons index','coords','euler char', 'boundary slopes', 'qtons surface', 'spinning slopes')
    treedata = [(i,verts[i].surface_index,verts[i].coords, verts[i].euler_char, str(verts[i].boundary_slopes), 'S_{},{}'.format(verts[i].genus, verts[i].num_boundary_comps), str(W.regina_bdy_slopes(verts[i].surface_index))) for i in vert_list]
    make_treeview_table(frame, treedata, column_names, with_scrollbar)

def make_facets_table(tnorm_app, frame, dim):
    if dim in [0,tnorm_app.ball.dimension]:
        make_vertices_table(tnorm_app, frame, 'all')
    else:
        W = tnorm_app.wrapper
        B = tnorm_app.ball
        facets = B.facets(dim)
        width = 7
        column_names = ('face','vertex','coords','euler char', 'boundary slopes', 'qtons surface')
        treedata = []
        for j in range(len(facets)):
            f = facets[j]
            vert_list = [v.index for v in f.vertices]
            verts = {i:B.vertices[i] for i in vert_list}
            v0 = vert_list[0]
            treedata.append((j,v0,verts[v0].coords, verts[v0].euler_char, str(verts[v0].boundary_slopes), 'S_{},{}'.format(verts[v0].genus, verts[v0].num_boundary_comps)))
            for i in vert_list[1:]:
                treedata.append(('',i,verts[i].coords, verts[i].euler_char, str(verts[i].boundary_slopes), 'S_{},{}'.format(verts[i].genus, verts[i].num_boundary_comps)))
            treedata.append(('','','','','','',''))
        make_treeview_table(frame, treedata, column_names)
    tnorm_app.stop_spin()

def make_dnb_vertices_table(tnorm_app, frame):
    W = tnorm_app.wrapper
    B = tnorm_app.ball
    DB = tnorm_app.dual_ball
    dim = B.dimension - 1
    facets = B.facets(dim)
    assert DB.num_vertices == len(facets)
    width = 7
    column_names = ('vertex', 'dual to face','vertices of dual face','coords','euler char', 'boundary slopes', 'qtons surface')
    treedata = []
    for j in range(len(facets)):
        f = facets[j]
        vert_list = [v.index for v in f.vertices]
        verts = {i:B.vertices[i] for i in vert_list}
        v0 = vert_list[0]
        treedata.append((j,j,v0,verts[v0].coords, verts[v0].euler_char, str(verts[v0].boundary_slopes), 'S_{},{}'.format(verts[v0].genus, verts[v0].num_boundary_comps)))
        for i in vert_list[1:]:
            treedata.append(('','',i,verts[i].coords, verts[i].euler_char, str(verts[i].boundary_slopes), 'S_{},{}'.format(verts[i].genus, verts[i].num_boundary_comps)))
        treedata.append(('','','','','','',''))
    make_treeview_table(frame, treedata, column_names)
    tnorm_app.stop_spin()


def make_treeview_table(frame, treedata, column_names, with_scrollbar=True):
    if with_scrollbar:
        xscrollbar = ttk.Scrollbar(frame, orient='horizontal')
        yscrollbar = ttk.Scrollbar(frame, orient='vertical')
        tree = ttk.Treeview(frame, columns=column_names, yscrollcommand=yscrollbar.set, xscrollcommand=xscrollbar.set, show=['headings'])
    else:
        tree = ttk.Treeview(frame, columns=column_names, show=['headings'])
    for i in range(len(column_names)):
        col = column_names[i]
        width = int(ceil(max([len(col), max([len(str(treedata[j][i])) for j in range(len(treedata))])])*7))
        tree.heading(col, text=col, anchor='center')
        tree.column(col, minwidth=0, width=width, stretch=tk.NO, anchor='center')
    for x in treedata:
        tree.insert('', 'end', values=x)
    frame.columnconfigure(0, weight=1) # column with treeview
    frame.rowconfigure(0, weight=1) # row with treeview  
    tree.grid(row=0, column=0, sticky='nswe')
    if with_scrollbar:
        xscrollbar.config(command=tree.xview)
        yscrollbar.config(command=tree.yview)
        yscrollbar.grid(row=0, column=1, sticky='ns')
        xscrollbar.grid(row=1, column=0, sticky='we')



