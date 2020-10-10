from __future__ import print_function
import tnorm

import tkinter as tk

from tkinter import ttk
from math import ceil



def make_qtons_table(tnorm_app):
#    try:
    qd = tnorm_app.wrapper.qtons_info
    treedata = [(i, qd[i]['image_in_H2'], qd[i]['euler_char'], 'S_{},{}'.format(qd[i]['genus'], qd[i]['num_boundary_comps']), qd[i]['is_embedded'], str(qd[i]['boundary_slopes'][0]), str(qd[i]['boundary_slopes'][1]), qd[i]['is_norm_minimizing'], qd[i]['over_facet']) for i in qd]
    column_names = ('qtons index', 'coords', 'euler char', 'topology', 'is_embedded', 'outward_bdy', 'inward_bdy','norm minimizing?','over facet')
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
        vert_list = range(len(B.vertices()))
    height = len(vert_list)
    verts = {i:B.vertices()[i] for i in vert_list}
    column_names = ('vertex', 'qtons index','coords','topology','euler char', 'is_embedded','outward_bdy','inward_bdy')
    treedata = [(i, verts[i].qtons_index(), verts[i].coords(), 'S_{},{}'.format(verts[i].genus(), verts[i].num_boundary_comps()), verts[i].euler_char(), verts[i].is_embedded(), str(verts[i].boundary_slopes()[0]), str(verts[i].boundary_slopes()[1])) for i in vert_list]
    make_treeview_table(frame, treedata, column_names, with_scrollbar)

def make_facets_table(tnorm_app, frame, dim):
    if dim in [0,tnorm_app.ball.dimension()]:
        make_vertices_table(tnorm_app, frame, 'all')
    else:
        W = tnorm_app.wrapper
        B = tnorm_app.ball
        facets = B.facets(dim)
        column_names = ('facet','vertex','qtons index','coords','topology','euler char','is_embedded', 'outward_bdy','inward_bdy')
        width = len(column_names)
        treedata = []
        for j in range(len(facets)):
            f = facets[j]
            vert_list = [v.index() for v in f.vertices()]
            verts = {i:B.vertices()[i] for i in vert_list}
            v0 = vert_list[0]
            treedata.append((j,v0, verts[v0].qtons_index(), verts[v0].coords(), 'S_{},{}'.format(verts[v0].genus(), verts[v0].num_boundary_comps()), verts[v0].euler_char(), verts[v0].is_embedded(), str(verts[v0].boundary_slopes()[0]), str(verts[v0].boundary_slopes()[1])))
            for i in vert_list[1:]:
                treedata.append(('',i, verts[i].qtons_index(), verts[i].coords(), 'S_{},{}'.format(verts[i].genus(), verts[i].num_boundary_comps()), verts[i].euler_char(), verts[i].is_embedded(), str(verts[i].boundary_slopes()[0]), str(verts[i].boundary_slopes()[1])))
            treedata.append(tuple('' for i in range(width)))
        make_treeview_table(frame, treedata, column_names)
    tnorm_app.stop_spin()

def make_dnb_vertices_table(tnorm_app, frame):
    W = tnorm_app.wrapper
    B = tnorm_app.ball
    DB = tnorm_app.dual_ball
    dim = B.dimension() - 1
    facets = B.facets(dim)
    assert DB.num_vertices() == len(facets)
    column_names = ('vertex', 'dual to face','vertices of dual face','qtons_index','coords', 'topology','euler char', 'outward_bdy', 'inward_bdy')
    width = len(column_names)
    treedata = []
    for j in range(len(facets)):
        f = facets[j]
        if dim == 0:
            vert_list = [f.index()]
        else:
            vert_list = [v.index() for v in f.vertices()]
        verts = {i:B.vertices()[i] for i in vert_list}
        v0 = vert_list[0]
        treedata.append((j,j,v0, verts[v0].qtons_index(), verts[v0].coords(), 'S_{},{}'.format(verts[v0].genus(), verts[v0].num_boundary_comps()), verts[v0].euler_char(), str(verts[v0].boundary_slopes())[0], str(verts[v0].boundary_slopes())[1]))
        for i in vert_list[1:]:
            treedata.append(('','',i, verts[i].qtons_index(), verts[i].coords(), 'S_{},{}'.format(verts[i].genus(), verts[i].num_boundary_comps()), verts[i].euler_char(), str(verts[i].boundary_slopes())[0], str(verts[v0].boundary_slopes())[1]))
        treedata.append(tuple('' for i in range(width)))
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



