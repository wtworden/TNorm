from __future__ import print_function

import tkinter as tk
from tkinter import ttk


def generate_summary(tnorm_app):
    try:
        print('Generating summary page... ', end='')
        W=tnorm_app.wrapper
        B=tnorm_app.ball
        P=tnorm_app.ball.polyhedron()
        M=tnorm_app.wrapper.manifold()
        try:
            sym_group = M.symmetry_group()
        except ValueError:
            sym_group = 'Failed'
        sections = ['Manifold', 'Thurston norm ball', 'Normal surfaces']
        names = [('rank H_2(M;bdy M)', 'volume', 'num tetrahedra', 'num cusps', 'symmetry group'),('num vertices', 'vertex degrees', 'vertices qtons rep (genus,punctures)', 'num faces', 'top dim faces (<verts>)'),('allows non-admissible', 'num qtons', 'num non-trivial in H2', 'num norm minimizing')]
        qtons_inds = [v.qtons_index() for v in B.vertices()]
        values = [(M.homology().rank(), M.volume(), M.num_tetrahedra(), M.num_cusps(), sym_group), (len(P.vertices()), [sum(P.vertex_adjacency_matrix()[i]) for i in range(P.n_vertices())], [(v.genus(),v.num_boundary_comps()) for v in B.vertices()], len(P.faces(P.dim()-1)), ['<{}>'.format(' '.join([str(v.index()) for v in face.vertices()])) for face in P.faces(P.dim()-1)]), (W.allows_non_admissible(), W.qtons().size(), len([1 for i in range(W.qtons().size()) if W.over_facet(i)!=None]), len([1 for i in range(W.qtons().size()) if W.is_norm_minimizing(i)]))]
        summary_page(tnorm_app.SummaryTab, sections, names, values)
        print('Done.')
    except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()

def generate_nb_overview(tnorm_app):
    try:
        W=tnorm_app.wrapper
        B=tnorm_app.ball
        P=tnorm_app.ball.polyhedron()
        M=tnorm_app.wrapper.manifold()

        sections = ['Overview']
        names = [('dimension', 'is compact?', 'num vertices', 'num rays', 'vertex degrees', 'vertices qtons rep (genus,punctures)', 'all vertices admissible?', 'num faces', 'faces (<vertices by index>)')]
        values = [(B.dimension(), B.is_compact(), len(P.vertices()), B.num_rays(), [sum(P.vertex_adjacency_matrix()[i]) for i in range(P.n_vertices())], [(v.genus(),v.num_boundary_comps()) for v in B.vertices()], B.all_vertices_admissible(), len(P.faces(P.dim()-1)), ['<{}>'.format(' '.join([str(B.index_of_poly_vert(v)) for v in face.vertices()])) for face in P.faces(P.dim()-1)])]
        summary_page(tnorm_app.NormBallTab, sections, names, values)
    except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()

def generate_dnb_overview(tnorm_app):
    try:
        W=tnorm_app.wrapper
        B=tnorm_app.dual_ball
        P=tnorm_app.dual_ball.polyhedron()
        M=tnorm_app.wrapper.manifold()

        sections = ['Overview']
        names = [('dimension', 'is full dimensional?', 'num vertices', 'vertex degrees', 'num faces', 'faces (<vertices by index>)')]
        values = [(B.dimension(), B.is_full_dimensional(), len(P.vertices()), [sum(P.vertex_adjacency_matrix()[i]) for i in range(P.n_vertices())], len(P.faces(P.dim()-1)), ['<{}>'.format(' '.join([str(B.index_of_poly_vert(v)) for v in face.vertices()])) for face in P.faces(P.dim()-1)])]
        summary_page(tnorm_app.DualNormBallTab, sections, names, values)
    except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()

def summary_page(frame, sections, names, values):

    col1_width = max([len('{}'.format(section_names[i])) for section_names in names for i in range(len(section_names))]) + 5
    col2_width = max([len('{}'.format(section_values[i])) for section_values in values for i in range(len(section_values))]) + 5
    
    text = tk.Text(frame, wrap='none', bg="#EAEAEA", height=0, borderwidth=-2, padx=-1, pady=-1, insertwidth=0)
    xscrollbar = ttk.Scrollbar(frame, orient='horizontal')
    yscrollbar = ttk.Scrollbar(frame, orient='vertical')
    xscrollbar.config(command=text.xview)
    yscrollbar.config(command=text.yview)
    text['xscrollcommand'] = xscrollbar.set
    text['yscrollcommand'] = yscrollbar.set
    frame.columnconfigure(0, weight=1) # column with treeview
    frame.rowconfigure(0, weight=1) # row with treeview  
    yscrollbar.grid(row=0, column=1, sticky='ns')
    text.grid(row=0, column=0, sticky='nswe')
    xscrollbar.grid(row=1, column=0, sticky='we')
    column1 = tk.Text(text, state='normal', wrap='none', width = 44, bg="#EAEAEA", borderwidth=-2, padx=-1, pady=-1)
    column2 = tk.Text(text, state='normal', wrap='none', width = col2_width, bg="#EAEAEA", borderwidth=-2, padx=-1, pady=-1)
    text.window_create(tk.END, window=column1)
    text.window_create(tk.END, window=column2)
    
    column1.tag_add('left', tk.INSERT)
    column2.tag_add('left', tk.INSERT)
    column1.tag_add('right', tk.INSERT)
    column1.tag_add('right', tk.INSERT)
    column1.tag_configure('left', justify='left')
    column2.tag_configure('left', justify='left')
    column1.tag_configure('right', justify='right')
    column1.tag_configure('right', justify='right')
    column1.bind("<1>", lambda event: column1.focus_set())
    column2.bind("<1>", lambda event: column2.focus_set())

    height = 1
    for i in range(len(sections)):
        column1.insert(tk.END, '\n\n   {}:{}\n'.format(sections[i], ' '*(col1_width-5-len(sections[i]))), ('left'))
        column2.insert(tk.END, '\n\n{}\n'.format(' '*(col2_width-2)), ('left'))
        height += 3
        for j in range(len(names[i])):
            height += 1
            column1.insert(tk.END, '{}:\n'.format(names[i][j]), ('right'))
            column2.insert(tk.END, ' {}\n'.format(values[i][j]), ('left'))
    column1.config(state='disable', height=height)
    column2.config(state='disable', height=height)






