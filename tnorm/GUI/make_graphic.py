from __future__ import print_function
import tnorm

import tkinter as tk
from tnorm.utilities import make_temp_directory
import os

from tkinter import ttk
import tempfile
from tnorm.GUI.canvas_image import CanvasImage


def make_canvas_image(tnorm_app, tab, filename):
    tnorm_app.ZoomFrame = tk.Frame(tab)
    tnorm_app.ZoomFrame.pack(side=tk.LEFT,expand=False, fill=tk.Y)
    tnorm_app.ImageFrame = tk.Frame(tab)  # placeholder of the ImageFrame object
    tnorm_app.ImageFrame.pack(side=tk.RIGHT, expand=True, fill=tk.BOTH)

    tnorm_app.ImageFrame.rowconfigure(0, weight=1)  # make the CanvasImage widget expandable
    tnorm_app.ImageFrame.columnconfigure(0, weight=1)
    canvas = CanvasImage(tnorm_app, filename)
    canvas.grid(row=0, column=0)


def make_hasse(tnorm_app, tab):
#    try:
    if tab == tnorm_app.NormBallTab:
        G = tnorm_app.ball.get_hasse_diagram
    elif tab == tnorm_app.DualNormBallTab:
        G = tnorm_app.dual_ball.get_hasse_diagram
    with make_temp_directory() as temp_dir:
        filename = os.path.join(temp_dir, 'hasse.png')
        G.save_image(filename)
        make_canvas_image(tnorm_app, tab, filename)

#    except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()


def _show_polyhedron(tnorm_app, online, tab):
    #try:
    if tab == tnorm_app.NormBallTab:
        ball = tnorm_app.ball
        dual = False
    elif tab == tnorm_app.DualNormBallTab:
        ball = tnorm_app.dual_ball
        dual = True
    if not ball.is_compact():
        P = ball.polyhedron_mod_rays()
    else:
        P = ball.polyhedron()
    if P.dim() in [3,4]:
        ball.plot(viewer='x3d', online=online)
    elif P.dim() in [1,2]:
        with make_temp_directory() as temp_dir:
            filename = os.path.join(temp_dir,'plot2d.png')
            plt = ball._plot2d(dual)
            plt.save(filename)
            make_canvas_image(tnorm_app, tab, filename)
    else:
        print('Error: can\'t plot a polyhedron of dimension > 4.')
    #except Exception as e: print('Error: {}'.format(e))
    tnorm_app.stop_spin()











