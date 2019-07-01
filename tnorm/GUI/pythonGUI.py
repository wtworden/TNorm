#!/usr/bin/python

import snappy
import threading
import Queue
import time

from tnorm.TN_wrapper import TN_wrapper
from tnorm.GUI.console import ConsoleText
from tnorm.hasse import get_hasse
from tnorm.GUI.zoom import *
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

class MyApp:
    def __init__(self, parent):

         #------ constants for controlling layout of buttons ------
        button_width = 12

        # -------------- end constants ----------------


        # -------------- end constants ----------------

        self.parent = parent
        self.parent.geometry("900x600")

        ### Our topmost frame is called myContainer1
        self.main_frame = tk.Frame(parent) ###
        self.main_frame.pack(expand=YES, fill=BOTH)


        ### We will use HORIZONTAL (left/right) orientation inside myContainer1.
        ### Inside myContainer1, we create control_frame and demo_frame.

        # header frame ------------------------------------------

        self.input_frame = tk.Frame(self.main_frame,bg="#EAEAEA",borderwidth=3)
        self.input_frame.pack(side=TOP, expand=NO, fill=X, padx=2, pady=2)

        # header frame: input box -------------------------------

        self.input_box_frame = tk.Frame(self.input_frame, bg="#EAEAEA")
        self.input_box_frame.pack(side=LEFT)

        tk.Label(self.input_box_frame, text='Manifold', bg="#EAEAEA").grid(row=0) 
        self.mfld_entry = tk.Entry(self.input_box_frame) 
        self.mfld_entry.bind("<Return>", self.load_on_enter)
        self.mfld_entry.grid(row=0, column=1)

        # header frame: load button -----------------------------

        self.load_button_frame = tk.Frame(self.input_frame, bg="#EAEAEA")
        self.load_button_frame.pack(side=LEFT)

        self.loadButton = ttk.Button(self.load_button_frame, text="Load", command=self.load)
        self.loadButton.pack(side=LEFT)

        # header frame: progress bar ----------------------------

        self.progress_bar_frame = tk.Frame(self.input_frame, bg="#EAEAEA")
        self.progress_bar_frame.pack(side=LEFT)



        # body frame --------------------------------------------
       
        self.bottom_frame = tk.Frame(self.main_frame,bg="#EAEAEA",borderwidth=3)
        self.bottom_frame.pack(side=BOTTOM, expand=YES, fill=BOTH,padx=2, pady=2)

        # body frame: control frame -----------------------------

        self.control_frame = tk.Frame(self.bottom_frame, borderwidth=3,  relief=FLAT, bg="#EAEAEA", width=200) ###
        self.control_frame.pack(side=LEFT, expand=NO, fill=Y, anchor=W)

        myMessage="Options"
        tk.Label(self.control_frame, text=myMessage, justify=LEFT, bg="#EAEAEA").pack(side=TOP, anchor=W)

        # body frame: control frame: buttons frame --------------

        self.buttons_frame = tk.Frame(self.control_frame, bg="#EAEAEA") ###
        self.buttons_frame.pack(side=TOP, expand=YES, fill=Y)

        self.vertexListButton = ttk.Button(self.buttons_frame, text="vertex list",command=self.showVertices,width=button_width)
        self.vertexListButton.pack()

        self.faceListButton = ttk.Button(self.buttons_frame, text="(n-1)-dim faces",width=button_width, command=self.showFaces)
        self.faceListButton.pack()

        self.facetButtonFrame = tk.Frame(self.buttons_frame, bg="#EAEAEA")
        self.facetButtonFrame.pack(side=TOP, expand=NO)

        self.dimEntry = tk.Entry(self.facetButtonFrame, width=2)
        self.dimEntry.pack(side=LEFT)

#        self.dimEntryDash = tk.Label(self.facetButtonFrame, text='-',bg="#EAEAEA")
#        self.dimEntryDash.pack(side=LEFT)

        self.facetButton = ttk.Button(self.facetButtonFrame, text='-dim facets', width=button_width-4, command=self.showFacets)
        self.facetButton.pack(side=LEFT)

        self.showPolyhedronButton = ttk.Button(self.buttons_frame, text="view norm ball",width=button_width, command=self.showPolyhedron)
        self.showPolyhedronButton.pack()

        self.showHasseButton = ttk.Button(self.buttons_frame, text="show hasse",width=button_width, command=self.showHasse)
        self.showHasseButton.pack()

        # body frame: control frame: cancel button frame --------------

        self.cancel_button_frame = tk.Frame(self.control_frame)
        self.cancel_button_frame.pack(side=BOTTOM, expand=NO, anchor=SW)

        self.cancelButton = ttk.Button(self.cancel_button_frame,
            text="Cancel",
            width=button_width,
            command=self.cancelButtonClick
            )
        self.cancelButton.pack(side=BOTTOM, anchor=S)

        # body frame: output frame ------------------------------

        self.out_frame = tk.Frame(self.bottom_frame, bg="#EAEAEA",borderwidth=3) ###
        self.out_frame.pack(side=RIGHT, expand=YES, fill=BOTH)

        # body frame: output frame: canvas frame ------------

        self.canvas_frame = tk.Frame(self.out_frame,bg="#EAEAEA",borderwidth=3,relief=RIDGE)
        self.canvas_frame.pack(side=TOP, expand=YES, fill=BOTH)  ###

        # body frame: output frame: console frame ----------

        self.console_frame = tk.Frame(self.out_frame,borderwidth=3,
            height=100, relief=RIDGE) ###
        self.console_frame.pack(side=TOP, fill=X)

        self.console = ConsoleText(self.console_frame, height=7, width=20)
        self.console.pack(expand=YES,fill=X)
        self.console.start()

    def clear_canvas(self):
        for widget in self.canvas_frame.winfo_children():
            widget.destroy()

    def showHasse(self):
        G = get_hasse(self.ball.polyhedron)
        temp_directory = tempfile.gettempdir()
        G.save_image(temp_directory+'/hasse.png')
        #image = PIL.Image.open(temp_directory+'/hasse.png')
        #self.img = ImageTk.PhotoImage(image)
        self.clear_canvas()

        self.zoom_frame = ttk.Frame(self.canvas_frame)
        self.zoom_frame.pack(side=tk.LEFT,expand=False, fill=tk.Y)
        self.image_frame = ttk.Frame(self.canvas_frame)  # placeholder of the ImageFrame object
        self.image_frame.pack(side=tk.RIGHT, expand=True, fill=tk.BOTH)

        self.image_frame.rowconfigure(0, weight=1)  # make the CanvasImage widget expandable
        self.image_frame.columnconfigure(0, weight=1)
        canvas = CanvasImage(self, temp_directory+'/hasse.png')
        canvas.grid(row=0, column=0)
#        canvas.pack(expand=YES,fill=BOTH)
#        canvas.create_image(0, 0, anchor=NW, image=self.img)

    def showPolyhedron(self):
        self.ball.plot()

    def showFaces(self):
        pass

    def showFacets(self):
        dim = self.dimEntry.get()
        P = self.ball.polyhedron

        pass

    def progress(self):
        self.prog_bar = ttk.Progressbar(
            self.progress_bar_frame, orient="horizontal",
            length=200, mode="indeterminate"
            )
        self.prog_bar.pack(side=TOP)

    def load_on_enter(self,event):
        self.load()

    def cancelButtonClick(self):
        self.parent.destroy()

    def showVertices(self):
        self.clear_canvas()
        if self.ball:
            W = self.wrapper
            B = self.ball
            height = B.numVertices
            width = 7
            verts = [B.verticesList.lookup(i) for i in range(height)]
            header = ['vertex','coords','euler char', 'H_1 boundary', 'qtons bdy', 'qtons surface','rep surface']
            for j in range(width): 
                b = Label(self.canvas_frame, text=header[j], bg="#EAEAEA")
                b.grid(row=0, column=j, ipadx=4)        
            for i in range(height): 
                v = verts[i]
                b = Label(self.canvas_frame, text='{}'.format(i), bg="#EAEAEA").grid(row=i+1, column=0)
                b = Label(self.canvas_frame, text='{}'.format(v.coords), bg="#EAEAEA").grid(row=i+1, column=1, ipadx=6)
                b = Label(self.canvas_frame, text='{}'.format(v.eulerChar), bg="#EAEAEA").grid(row=i+1, column=2)
                b = Label(self.canvas_frame, text='{}'.format(v.H1BoundarySlopes), bg="#EAEAEA").grid(row=i+1, column=3,ipadx=6)
                b = Label(self.canvas_frame, text='{}'.format(v.spinningSlopes), bg="#EAEAEA").grid(row=i+1, column=4,ipadx=6)
                b = Label(self.canvas_frame, text='{}'.format('---'), bg="#EAEAEA").grid(row=i+1, column=5)
                b = Label(self.canvas_frame, text='S_{},{}'.format(v.genus,v.numPunctures), bg="#EAEAEA").grid(row=i+1, column=6)

        else:
            print('please load a manifold first.')    


    def compute_TN(self,string):
        try:
            M=snappy.Manifold(string)
            result = {}
            t = threading.Thread(target=self.TN_wrapper,args=(M,result))
            t.start()
        except:
            print("invalid entry")

    def TN_wrapper(self, M, result):
        result[0] = TN_wrapper(M)
        print("Computing Thurston norm ball.\n")
        result[1] = result[0].normBall()
        self.wrapper = result[0]
        self.ball = result[1]
        self.showVertices()

    def load(self):
        string = self.mfld_entry.get() 
        self.compute_TN(string)

class ResizingCanvas(Canvas):
    def __init__(self,parent,**kwargs):
        Canvas.__init__(self,parent,**kwargs)
        self.bind("<Configure>", self.on_resize)
        self.height = self.winfo_reqheight()
        self.width = self.winfo_reqwidth()

    def on_resize(self,event):
        # determine the ratio of old width/height to new width/height
        wscale = float(event.width)/self.width
        hscale = float(event.height)/self.height
        self.width = event.width
        self.height = event.height
        # resize the canvas 
        self.config(width=self.width, height=self.height)
        # rescale all the objects tagged with the "all" tag
        self.scale("all",0,0,wscale,hscale)


root = tk.Tk()
root.title('tnorm')
myapp = MyApp(root)
root.mainloop()






