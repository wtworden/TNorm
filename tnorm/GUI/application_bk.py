#!/usr/bin/python

import snappy
import threading
import Queue
import time
import webbrowser


from tnorm.GUI.make_table import make_qtons_table, make_vertices_table
from tnorm.GUI.make_graphic import make_hasse, _show_polyhedron
from tnorm.utilities import *
from tnorm.TN_wrapper import TN_wrapper
from tnorm.GUI.console import ConsoleText
from tnorm.hasse import get_hasse
from tnorm.GUI.canvas_image import CanvasImage
from tnorm.GUI.summary_page import generate_summary
# background="..." doesn't work...

import PIL
from PIL import ImageTk

import tempfile
import re
import os
import io
import sys
#from math import sin, cos, pi, ceil, sqrt

FileType = io.TextIOWrapper if sys.version_info >= (3, 0) else file

try:
    import Tkinter as tk
    import tkFileDialog
    import tkMessageBox
except ImportError:  # Python 3.
    try:
        import tkinter as tk
        import tkinter.filedialog as tkFileDialog
        import tkinter.messagebox as tkMessageBox
    except ImportError:
        raise ImportError('Tkinter not available.')

try:
    import ttk as ttk
except ImportError:  # Python 3.
    try:
        from tkinter import ttk
    except ImportError:
        raise ImportError('Ttk not available.')

from tkinter import *

# Some constants.
if sys.platform in ['darwin']:
    COMMAND = {
        'new': 'Command+N',
        'open': 'Command+O',
        'save': 'Command+S',
        'close': 'Command+W',
        'lamination': 'Command+L',
        'erase': 'Command+D',
        'twist': 'Command+T',
        'halftwist': 'Command+H',
        'isometry': 'Command+I',
        'compose': 'Command+M'
        }
    COMMAND_KEY = {
        'new': '<Command-n>',
        'open': '<Command-o>',
        'save': '<Command-s>',
        'close': '<Command-w>',
        'lamination': '<Command-l>',
        'erase': '<Command-d>',
        'twist': '<Command-t>',
        'halftwist': '<Command-h>',
        'isometry': '<Command-i>',
        'compose': '<Command-m>'
        }
else:
    COMMAND = {
        'new': 'Ctrl+N',
        'open': 'Ctrl+O',
        'save': 'Ctrl+S',
        'close': 'Ctrl+W',
        'lamination': 'Ctrl+L',
        'erase': 'Ctrl+D',
        'twist': 'Ctrl+T',
        'halftwist': 'Ctrl+H',
        'isometry': 'Ctrl+I',
        'compose': 'Ctrl+M'
        }
    COMMAND_KEY = {
        'new': '<Control-n>',
        'open': '<Control-o>',
        'save': '<Control-s>',
        'close': '<Control-w>',
        'lamination': '<Control-l>',
        'erase': '<Control-d>',
        'twist': '<Control-t>',
        'halftwist': '<Control-h>',
        'isometry': '<Control-i>',
        'compose': '<Control-m>'
        }

class TNormApp:
    def __init__(self, parent):

         #------ constants for controlling layout of buttons ------
        BUTTON_WIDTH = 12
        BUTTON_YPAD = .1


        self.parent = parent

        ### Our topmost frame is called MainFrame
        self.MainFrame = tk.Frame(parent) ###
        self.MainFrame.pack(expand=YES, fill=BOTH)

        #  control frame -----------------------------

        self.ControlFrame = tk.Frame(self.MainFrame, borderwidth=3,  relief=FLAT, bg="#EAEAEA", width=200) ###
        self.ControlFrame.pack(side=LEFT, expand=NO, fill=Y, anchor=W)

        

        # IO frame ------------------------------

        self.IOFrame = tk.Frame(self.MainFrame, bg="#EAEAEA",borderwidth=3) ###
        self.IOFrame.pack(side=RIGHT, expand=YES, fill=BOTH)

        # IO frame: header frame ------------------------------------------

        self.InputFrame = tk.Frame(self.IOFrame,bg="#EAEAEA",borderwidth=3)
        self.InputFrame.pack(side=TOP, expand=NO, fill=X, padx=2, pady=2)

        # IO frame: header frame: input box -------------------------------

        self.InputBoxFrame = tk.Frame(self.InputFrame, bg="#EAEAEA")
        self.InputBoxFrame.pack(side=LEFT)

        tk.Label(self.InputBoxFrame, text='Manifold', bg="#EAEAEA").grid(row=0) 
        self.ManifoldEntry = tk.Entry(self.InputBoxFrame) 
        self.ManifoldEntry.bind("<Return>", self.load_on_enter)
        self.ManifoldEntry.grid(row=0, column=1)

        # IO frame: header frame: load button -----------------------------

        self.LoadButtonFrame = tk.Frame(self.InputFrame, bg="#EAEAEA")
        self.LoadButtonFrame.pack(side=LEFT)

        self.LoadButton = ttk.Button(self.LoadButtonFrame, text="Load", command=self.load)
        self.LoadButton.pack(side=LEFT)

        # control frame : help frame --------------------------------------

        self.HelpFrame = tk.Frame(self.ControlFrame, bg="#EAEAEA",height=30)
        self.HelpFrame.pack(side=TOP)

#        self.HelpButton


        # control frame: buttons frame -------------------------------------

        self.ButtonsFrame = tk.Frame(self.ControlFrame, bg="#EAEAEA") ###
        self.ButtonsFrame.pack(side=TOP, expand=YES, fill=Y)

        tk.Label(self.ButtonsFrame, text="Options", justify=LEFT, bg="#EAEAEA").pack(side=TOP, anchor=W)

        self.SummaryButton = ttk.Button(self.ButtonsFrame, text="summary",command=self.show_summary,width=BUTTON_WIDTH)
        self.SummaryButton.pack(pady = BUTTON_YPAD)

        self.spacer1 = tk.Canvas(self.ButtonsFrame, height=10, bg="#EAEAEA", width=100, highlightthickness=-2)
        self.spacer1.pack()

        self.VertexListButton = ttk.Button(self.ButtonsFrame, text="vertex list",command=self.show_vertices,width=BUTTON_WIDTH)
        self.VertexListButton.pack(pady = BUTTON_YPAD)

        self.FaceListButton = ttk.Button(self.ButtonsFrame, text="(n-1)-dim faces",width=BUTTON_WIDTH, command=self.show_faces)
        self.FaceListButton.pack(pady = BUTTON_YPAD)

        self.FacetButtonFrame = tk.Frame(self.ButtonsFrame, bg="#EAEAEA")
        self.FacetButtonFrame.pack(side=TOP, expand=NO)

        self.DimEntry = tk.Entry(self.FacetButtonFrame, width=2)
        self.DimEntry.pack(side=LEFT)

        self.FacetButton = ttk.Button(self.FacetButtonFrame, text='-dim facets', width=BUTTON_WIDTH-4, command=self.show_facets)
        self.FacetButton.pack(side=LEFT, pady = BUTTON_YPAD)

        self.spacer2 = tk.Canvas(self.ButtonsFrame, height=10, bg="#EAEAEA", width=100, highlightthickness=-2)
        self.spacer2.pack()

        self.ShowPolyhedronButton = ttk.Button(self.ButtonsFrame, text="view norm ball",width=BUTTON_WIDTH, command=self.show_polyhedron)
        self.ShowPolyhedronButton.pack(pady = BUTTON_YPAD)

        self.ShowHasseButton = ttk.Button(self.ButtonsFrame, text="show hasse diag.",width=BUTTON_WIDTH, command=self.show_hasse)
        self.ShowHasseButton.pack(pady = BUTTON_YPAD)

        self.QtonsInfo = ttk.Button(self.ButtonsFrame, text="qtons info",width=BUTTON_WIDTH, command=self.show_all_qtons)
        self.QtonsInfo.pack(pady = BUTTON_YPAD)


        # control frame: cancel button frame --------------

        self.CancelButtonFrame = tk.Frame(self.ControlFrame, bg="#EAEAEA")
        self.CancelButtonFrame.pack(side=BOTTOM, expand=NO, anchor=SW)

        self.CancelButton = ttk.Button(self.CancelButtonFrame, text="Cancel", width=6, command=self.cancel_button_click)
        self.CancelButton.pack(side=BOTTOM, pady=10)


        # output frame --------------------------------------------
       
        self.OutputFrame = tk.Frame(self.IOFrame,bg="#EAEAEA",borderwidth=3)
        self.OutputFrame.pack(side=BOTTOM, expand=YES, fill=BOTH,padx=2, pady=2)


        # IO frame: output frame: canvas frame --------------------

        self.CanvasFrame = tk.Frame(self.OutputFrame,bg="#EAEAEA",borderwidth=3,relief=RIDGE)
        self.CanvasFrame.pack(side=TOP, expand=YES, fill=BOTH)  ###

        # body frame: output frame: console frame -----------------

        self.StatusFrame = tk.Frame(self.OutputFrame,borderwidth=0,
            height=120, bg="#EAEAEA") ###
        self.StatusFrame.pack(side=TOP, fill=X)

        #  console frame: icon frame -------------------------------

        self.IconFrame = tk.Frame(self.StatusFrame, bg="#EAEAEA")
        self.IconFrame.pack(side=RIGHT)

        datadir = os.path.dirname(__file__)
        path = os.path.join(datadir, 'images', 'spinpoly', 'poly{}.gif')
        self.icon_frames = [PhotoImage(file=path.format(i)) for i in range(9)]

        self.IconGif = Label(self.IconFrame,image=self.icon_frames[0], bd=-2)
        self.IconGif.pack()

        # console ---------------------------------------------------

        self.ConsoleFrame = tk.Frame(self.StatusFrame, borderwidth=2,
            height=120, relief=RIDGE) ###
        self.ConsoleFrame.pack(side=LEFT, expand=YES, fill=X)

        self.console = ConsoleText(self.ConsoleFrame, height=8, width=20)
        self.console.pack(side=LEFT, expand=YES,fill=X)
        self.console.start()
        self.ball = None




    # clear window to load new table/graphic/etc ------------------------

    def clear_canvas(self):
        for widget in self.CanvasFrame.winfo_children():
            widget.destroy()


    # spinning polyhedron -----------------------------------------------

    def start_spin(self):
        self.computing = True
        self.IconFrame.after(0, self.update, 0)
        
    def update(self,ind):
        if self.computing == True:
            frame = self.icon_frames[ind]
            ind = (ind+1)%9
            self.IconGif.configure(image=frame)
            self.IconFrame.after(100, self.update, ind)
        else:
            return None

    def stop_spin(self):
        self.computing = False


    # load button commands -------------------------------------------

#    def load(self):
#        string = self.ManifoldEntry.get() 
#        self.compute_TN(string)


    def load(self):
        string = self.ManifoldEntry.get()
        try:
            M=snappy.Manifold(string)
            computing=True
            t = threading.Thread(target=self.load_computations,args=(M,))
            t.start()
            self.start_spin()

        except IOError:
            print("invalid entry")

    def load_computations(self, M):
        self.wrapper = TN_wrapper(M)
        self.ball = self.wrapper.TNormBall
        generate_summary(self)
        self.stop_spin()


    # option button commands ----------------------------------------------------

    def show_summary(self):
        computing=True
        self.clear_canvas()
        t = threading.Thread(target=generate_summary,args=(self,))
        t.start()
        self.start_spin()

    def show_vertices(self):
        if self.ball:
            self.clear_canvas()
            computing = True
            t = threading.Thread(target=make_vertices_table, args=(self,))
            t.start()
            self.start_spin()

        else:
            print('please load a manifold first.')    

    def show_faces(self):
        P = self.ball.polyhedron
        dim = P.dim()-1
        faces = P.faces(dim)


    def show_facets(self):
        dim = self.DimEntry.get()
        P = self.ball.polyhedron
        faces = P.faces(dim)
        column_names = ('face', 'vertex coords', )

        pass


    def show_all_qtons(self):
        computing = True
        self.clear_canvas()
        t = threading.Thread(target=make_qtons_table,args=(self,))
        t.start()
        self.start_spin()

    def show_hasse(self):
        computing = True
        self.clear_canvas()
        t = threading.Thread(target=make_hasse, args=(self,))
        t.start()
        self.start_spin()

    def show_schlegel(self):
        computing = True
        self.clear_canvas()
        t = threading.Thread(target=make_schlegel, args=(self,))
        t.start()
        self.start_spin()

    def show_polyhedron(self):
        computing = True
        online=False
        t = threading.Thread(target=_show_polyhedron, args=(self, online))
        t.start()
        self.start_spin()


    def cancel_button_click(self):
        self.parent.destroy()


    # button bindings for Enter key ----------------------------------

    def load_on_enter(self,event):
        self.load()


def start():
    root = tk.Tk()
    root.minsize(500,500)
    root.geometry("900x600")
    root.title('tnorm')
    myapp = TNormApp(root)
    datadir = os.path.dirname(__file__)
    icon_path = os.path.join(datadir, 'images', 'icon', 'icon.gif')
    img = tk.PhotoImage(file=icon_path)
    root.tk.call('wm', 'iconphoto', root._w, img)
    root.mainloop()



if __name__ == '__main__':
    start()

