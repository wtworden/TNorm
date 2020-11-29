#!/usr/bin/python

import threading
import platform
import os
import io
import sys

from tnorm.GUI.make_graphic import make_hasse, _show_polyhedron
from tnorm.GUI.make_table import make_qtons_table, make_all_vertices_table, make_facets_table, make_dnb_vertices_table
from tnorm.TN_wrapper import TN_wrapper
from tnorm.GUI.console import ConsoleText
from tnorm.utilities.hasse import get_hasse
from tnorm.GUI.canvas_image import CanvasImage
from tnorm.GUI.summary_page import generate_summary, generate_nb_overview, generate_dnb_overview
# background="..." doesn't work...

#import PIL
#from PIL import ImageTk
import snappy


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
        }
    COMMAND_KEY = {
        'new': '<Command-n>',
        'open': '<Command-o>',
        'save': '<Command-s>',
        'close': '<Command-w>',
        }
else:
    COMMAND = {
        'new': 'Ctrl+N',
        'open': 'Ctrl+O',
        'save': 'Ctrl+S',
        'close': 'Ctrl+W',
        }
    COMMAND_KEY = {
        'new': '<Control-n>',
        'open': '<Control-o>',
        'save': '<Control-s>',
        'close': '<Control-w>',
        }

class TNormApp:
    def __init__(self, parent):

         #------ constants for controlling layout of buttons ------
        BUTTON_WIDTH = 12
        BUTTON_YPAD = 0
        BG_COLOR = "#EAEAEA"


        self.parent = parent

        ### Our topmost frame is called MainFrame
        self.MainFrame = tk.Frame(parent) ###
        self.MainFrame.pack(expand=YES, fill=BOTH)

        #  main frame: control frame -----------------------------

        self.ControlFrame = tk.Frame(self.MainFrame, borderwidth=3,  relief=FLAT, bg=BG_COLOR, width=200) ###
        self.ControlFrame.pack(side=LEFT, expand=NO, fill=Y, anchor=W)

        # main frame: output frame ------------------------------

        self.OutputFrame = tk.Frame(self.MainFrame, bg=BG_COLOR,borderwidth=0) ###
        self.OutputFrame.pack(side=RIGHT, expand=YES, fill=BOTH)


        # control frame: load frame ------------------------------------------

        self.LoadFrame = tk.Frame(self.ControlFrame,bg=BG_COLOR,borderwidth=3)
        self.LoadFrame.pack(side=TOP, expand=NO, fill=X, padx=2, pady=2)

        # control frame: load frame: input box frame -------------------------------

        self.InputBoxFrame = tk.Frame(self.LoadFrame, bg=BG_COLOR)
        self.InputBoxFrame.pack(side=TOP)

        self.ManifoldLabel = tk.Label(self.InputBoxFrame, text='Manifold:', bg=BG_COLOR)
        self.ManifoldLabel.grid(row=0,column=0,sticky='w') 
        self.ManifoldEntry = tk.Entry(self.InputBoxFrame, width=10) 
        self.ManifoldEntry.bind("<Return>", self.load_on_enter)
        self.ManifoldEntry.grid(row=1, column=0)

        self.LoadButton = ttk.Button(self.InputBoxFrame, text="Load", command=self.load, padding=-8)
        self.LoadButton.grid(row=1,column=1)

        # control frame: load frame: load options frame ----------------------

        self.LoadOptionsFrame = tk.Frame(self.LoadFrame, bg=BG_COLOR)
        self.LoadOptionsFrame.pack(side=TOP,fill=X)

        self.allow_na_var = IntVar(None, 0)
        self.AllowsAdmissible = tk.Checkbutton(self.LoadOptionsFrame, bg=BG_COLOR, text='Allow non-admissible ', variable=self.allow_na_var)
        self.AllowsAdmissible.pack(side=TOP)

        self.force_simplicial = IntVar(None, 0)
        self.ForceSimplicial = tk.Checkbutton(self.LoadOptionsFrame, bg=BG_COLOR, text='Use simplicial H2 map', variable=self.force_simplicial)
        self.ForceSimplicial.pack(side=TOP)

        #self.BasisOption = tk.Frame(self.LoadOptionsFrame, bg=BG_COLOR)
        #self.BasisOption.pack(side=TOP, anchor=NW)
#
        #self.BasisLabel = tk.Label(self.BasisOption, bg=BG_COLOR, text='Peripheral curves:')
        #self.BasisLabel.grid(row=0,column=0,sticky='nw')
#
        #self.bdy_H1_basis_var = StringVar(None,'natural')
        #self.NaturalRadio = tk.Radiobutton(self.BasisOption, background=BG_COLOR, text='canonical (if link)', variable=self.bdy_H1_basis_var, value='natural')
        #self.NaturalRadio.grid(row=1, column=0, sticky='w', padx=10)
#
        #self.ShortestRadio = tk.Radiobutton(self.BasisOption, background=BG_COLOR, text='shortest', variable=self.bdy_H1_basis_var, value='shortest')
        #self.ShortestRadio.grid(row=2, column=0, sticky='w', padx=10)
        #self.bdy_H1_basis_var.set('natural')




        # control frame: buttons frame -------------------------------------
        self.ButtonsFrame = tk.Frame(self.ControlFrame, bg=BG_COLOR, borderwidth=0) ###
        self.ButtonsFrame.pack(side=TOP, expand=YES, fill=Y)


#        # control frame: buttons frame: view frame ------------------------------
#        self.ViewFrame = tk.Frame(self.ButtonsFrame, bg=BG_COLOR)
#        self.ViewFrame.pack(side=TOP, expand=NO, fill=X)
#
#        tk.Label(self.ViewFrame, text="View:", justify=LEFT, bg=BG_COLOR).pack(side=TOP, anchor=W)
#
#        self.SummaryButton = ttk.Button(self.ViewFrame, text="summary",command=self.show_summary,width=BUTTON_WIDTH)
#        self.SummaryButton.pack(pady = BUTTON_YPAD)
#
#        self.NormBallButton = ttk.Button(self.ViewFrame, text="norm ball",width=BUTTON_WIDTH)
#        self.NormBallButton.pack(pady = BUTTON_YPAD)
#
#        self.DualNormBallButton = ttk.Button(self.ViewFrame, text="dual norm ball",width=BUTTON_WIDTH)
#        self.DualNormBallButton.pack(pady = BUTTON_YPAD)
#        
#        self.QtonsInfo = ttk.Button(self.ViewFrame, text="qtons",width=BUTTON_WIDTH, command=self.show_all_qtons)
#        self.QtonsInfo.pack(pady = BUTTON_YPAD)
        # -----------------------------------------------------------------------


        #spacer ---------------------------------------------------------------------
        self.spacer1 = tk.Canvas(self.ButtonsFrame, height=10, bg=BG_COLOR, width=100, highlightthickness=-2)
        self.spacer1.pack()
        #----------------------------------------------------------------------------

        #control frame: buttons frame: options frame --------------------------------

        self.OptionsFrame = tk.Frame(self.ButtonsFrame, bg=BG_COLOR)
        self.OptionsFrame.pack(side=TOP, expand=YES, fill=BOTH)

        tk.Label(self.OptionsFrame, text="Options:", justify=LEFT, bg=BG_COLOR).pack(side=TOP, anchor=W)

        self.OverviewButton = ttk.Button(self.OptionsFrame, text="overview",command=self.show_overview,width=BUTTON_WIDTH)
        self.OverviewButton.pack(pady = BUTTON_YPAD)        

        self.VertexListButton = ttk.Button(self.OptionsFrame, text="vertex list",command=self.show_vertices,width=BUTTON_WIDTH)
        self.VertexListButton.pack(pady = BUTTON_YPAD)

        self.FaceListButton = ttk.Button(self.OptionsFrame, text="(n-1)-dim faces",width=BUTTON_WIDTH, command=self.show_faces)
        self.FaceListButton.pack(pady = BUTTON_YPAD)

        # control frame: buttons frame: options frame: facet button frame --------------
        self.FacetButtonFrame = tk.Frame(self.OptionsFrame, bg=BG_COLOR)
        self.FacetButtonFrame.pack(side=TOP, expand=NO)

        self.DimEntry = tk.Entry(self.FacetButtonFrame, width=2)
        self.DimEntry.pack(side=LEFT)

        self.FacetButton = ttk.Button(self.FacetButtonFrame, text='-dim facets', width=BUTTON_WIDTH-4, command=self.show_facets)
        self.FacetButton.pack(side=LEFT, pady = BUTTON_YPAD)
        # -------------------------------------------------------------------------------

#        self.spacer2 = tk.Canvas(self.ButtonsFrame, height=10, bg=BG_COLOR, width=100, highlightthickness=-2)
#        self.spacer2.pack()

        self.PolyhedronButton = ttk.Button(self.OptionsFrame, text="polyhedron",width=BUTTON_WIDTH, command=self.show_polyhedron)
        self.PolyhedronButton.pack(pady = BUTTON_YPAD)

        self.HasseButton = ttk.Button(self.OptionsFrame, text="hasse diagram",width=BUTTON_WIDTH, command=self.show_hasse)
        self.HasseButton.pack(pady = BUTTON_YPAD)




        # control frame: cancel button frame --------------

#        self.CancelButtonFrame = tk.Frame(self.ControlFrame, bg=BG_COLOR)
#        self.CancelButtonFrame.pack(side=BOTTOM, expand=NO, anchor=SW)

#        self.CancelButton = ttk.Button(self.CancelButtonFrame, text="Cancel", width=6, command=self.cancel_button_click)
#        self.CancelButton.pack(side=BOTTOM, pady=10)


        # main frame: output frame: notebook frame --------------------

        self.NotebookFrame = tk.Frame(self.OutputFrame, bg=BG_COLOR, pady=4)
        self.NotebookFrame.pack(side=TOP, expand=YES, fill=BOTH)  ###

        self.Notebook = Notebook = ttk.Notebook(self.NotebookFrame, padding=[0,0,7,0])

        self.SummaryTab = tk.Frame(self.NotebookFrame, bd=1, background='#DADADA')
        self.NormBallTab = tk.Frame(self.NotebookFrame, bd=1, background='#DADADA')
        self.DualNormBallTab = tk.Frame(self.NotebookFrame, bd=1, background='#DADADA')
        self.QtonsTab = tk.Frame(self.NotebookFrame, bd=1, background='#DADADA')
        self.HelpTab = tk.Frame(self.NotebookFrame, bd=1, background='#DADADA')
        Notebook.add(self.SummaryTab, text='Summary', sticky='nsew', padding=[0,0,0,0], state='disabled')
        Notebook.add(self.NormBallTab, text='Norm ball', sticky='nsew', padding=[0,0,0,0], state='disabled')
        Notebook.add(self.DualNormBallTab, text='Dual norm ball', sticky='nsew', padding=[0,0,0,0], state='disabled')
        Notebook.add(self.QtonsTab, text='qton surfaces', sticky='nsew', padding=[0,0,0,0], state='disabled')
        Notebook.add(self.HelpTab, text='Help', sticky='nsew', padding=[0,0,0,0], state='normal')
        Notebook.grid(row=0, column=0, sticky=tk.NSEW, padx=0, pady=0)
        self.NotebookFrame.grid_columnconfigure(0, weight=1)
        self.NotebookFrame.grid_rowconfigure(0, weight=1)
        self.disable_all_buttons()
        Notebook.bind('<<NotebookTabChanged>>', self.update_current_tab)
        self.tabs = [self.SummaryTab, self.NormBallTab, self.DualNormBallTab, self.QtonsTab, self.HelpTab]

        # main frame: output frame: status frame -----------------

        self.StatusFrame = tk.Frame(self.OutputFrame,borderwidth=0,
            height=120, bg=BG_COLOR, pady=6) ###
        self.StatusFrame.pack(side=TOP, fill=X)

        # main frame: output frame: buttons frame: icon frame -------------------------------

        self.IconFrame = tk.Frame(self.ButtonsFrame, bg=BG_COLOR)
        self.IconFrame.pack(side=TOP, expand=YES)

        datadir = os.path.dirname(__file__)
        path = os.path.join(datadir, 'images', 'spinpoly', 'poly{}.gif')
        self.icon_frames = [PhotoImage(file=path.format(i)) for i in range(9)]

        self.IconGif = Label(self.IconFrame,image=self.icon_frames[0], bd=-2)
        self.IconGif.pack(side=TOP)

        # main frame: output frame: status frame: console frame ---------------------------------------------------

        self.ConsoleFrame = tk.Frame(self.StatusFrame, borderwidth=2,
            height=120, relief=RIDGE) ###
        self.ConsoleFrame.grid(row=0,column=0, sticky='nsew')

        self.ConsolePadding = tk.Frame(self.StatusFrame, width=8, bg=BG_COLOR)
        self.ConsolePadding.grid(row=0,column=1)
        self.StatusFrame.grid_columnconfigure(0, weight=1)

        self.console = ConsoleText(self.ConsoleFrame, height=4, width=20)
        self.console.pack(side=LEFT, expand=YES,fill=X)
        self.console.start()
        self.ball = None
        self.dual_ball = None

        s = ttk.Style()
        s.configure("Treeview.Heading", relief='flat')

        if platform.system() == 'Darwin':
            s.theme_settings("default", {"TNotebook.Tab": {"configure": {"padding": [30, 30]}}})
        s.map("TButton",
            fieldbackground=[('disabled', '#F0F0F0')], foreground=[('disabled', '#909090')])
        s.map("TNotebook.Tab",
            fieldbackground=[('disabled', '#F0F0F0')], foreground=[('disabled', '#909090')])

        self.show_welcome()

    def current_tab(self):
        if self.Notebook.select() != '':
            tab_name = self.Notebook.tab(self.Notebook.select(),'text')
            if tab_name == 'Summary':
                return self.SummaryTab
            elif tab_name == 'Norm ball':
                return self.NormBallTab
            elif tab_name == 'Dual norm ball':
                return self.DualNormBallTab
            elif tab_name == 'qton surfaces':
                return self.QtonsTab
            elif tab_name == 'Help':
                return self.HelpTab
        else:
            return None

    def update_current_tab(self, event=None):
        self.parent.update_idletasks()
        tab = self.current_tab()
        if tab == self.SummaryTab:
            self.disable([self.DimEntry, self.OverviewButton, self.VertexListButton, self.FaceListButton, self.FacetButton, self.PolyhedronButton, self.HasseButton])            
        if tab == self.NormBallTab:
            self.DimEntry.delete(0,END)
            self.enable([self.DimEntry, self.OverviewButton, self.VertexListButton, self.FaceListButton, self.FacetButton, self.PolyhedronButton, self.HasseButton])
            if len(self.NormBallTab.winfo_children()) != 0:
                pass
            else:
                self.show_overview()
        elif tab == self.DualNormBallTab:
            self.enable([self.DimEntry, self.OverviewButton, self.VertexListButton, self.FaceListButton, self.FacetButton, self.PolyhedronButton, self.HasseButton])
            self.disable([self.DimEntry, self.FaceListButton, self.FacetButton])
            self.DimEntry.delete(0,END)
            if len(self.DualNormBallTab.winfo_children()) != 0:
                pass
            else:
                self.show_overview()

        elif tab == self.QtonsTab:
            self.disable([self.DimEntry, self.OverviewButton, self.VertexListButton, self.FaceListButton, self.FacetButton, self.PolyhedronButton, self.HasseButton])
            if len(self.QtonsTab.winfo_children()) == 0:
                self.show_all_qtons()
            else:
                pass
        elif tab == self.HelpTab:
            self.disable([self.DimEntry, self.VertexListButton, self.FaceListButton, self.FacetButton, self.PolyhedronButton, self.HasseButton])

        self.parent.update_idletasks()

    # disable tabs and option buttons

    def disable_all_tabs(self):
        current_tab = self.current_tab()
        # need to disable the current tab last, to avoid invoking update_current_tab().
        for tab in self.tabs:
            if tab != current_tab:
                self.Notebook.tab(tab, state='disabled')
        if current_tab != None:
            self.Notebook.tab(current_tab, state='disabled')


    def disable_all_buttons(self):
        self.VertexListButton.configure(state='disabled')
        self.FaceListButton.configure(state='disabled')
        self.FacetButton.configure(state='disabled')
        self.PolyhedronButton.configure(state='disabled')
        self.HasseButton.configure(state='disabled')
        self.DimEntry.configure(state='disabled')        

    def tabs_normal(self):
        self.Notebook.tab(self.SummaryTab, state='normal')
        self.Notebook.tab(self.NormBallTab, state='normal')
        self.Notebook.tab(self.DualNormBallTab, state='normal')
        self.Notebook.tab(self.QtonsTab, state='normal')
        self.Notebook.tab(self.HelpTab, state='normal')

    def enable(self, widgets):
        for widget in widgets:
            widget.configure(state='normal')
    
    def disable(self, widgets):
        for widget in widgets:
            widget.configure(state='disabled')


    # clear tab to load new table/graphic/etc ------------------------

    def clear_all_tabs(self):
        for tab in [self.SummaryTab, self.NormBallTab, self.DualNormBallTab, self.QtonsTab]:
            self.clear_tab(tab)

    def clear_tab(self, tab):
        for widget in tab.winfo_children():
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

    def load(self):
        try:
            M_str = self.ManifoldEntry.get()
            computing=True
            M = snappy.Manifold(M_str)
            t = threading.Thread(target=self.load_computations, args=(M,))
            t.start()
            self.start_spin()

        except IOError:
            print("invalid entry")

    def load_computations(self, M):
#        try:
        #self.disable([self.LoadButton, self.ManifoldEntry, self.NaturalRadio, self.ShortestRadio, self.AllowsAdmissible, self.BasisLabel, self.ForceSimplicial, self.OverviewButton])
        self.disable([self.LoadButton, self.ManifoldEntry, self.AllowsAdmissible, self.ForceSimplicial, self.OverviewButton])
        self.disable_all_buttons()
        self.disable_all_tabs()
        self.clear_all_tabs()
        self.wrapper = TN_wrapper(M, allows_non_admissible=bool(self.allow_na_var.get()), force_simplicial_homology=bool(self.force_simplicial.get()))
        self.ball = self.wrapper.norm_ball
        self.dual_ball = self.wrapper.dual_norm_ball
        generate_summary(self)
        self.tabs_normal()
        #self.enable([self.LoadButton, self.ManifoldEntry, self.NaturalRadio, self.ShortestRadio, self.AllowsAdmissible, self.ForceSimplicial, self.BasisLabel])
        self.enable([self.LoadButton, self.ManifoldEntry, self.AllowsAdmissible, self.ForceSimplicial])
        self.Notebook.select(self.SummaryTab)
#        except Exception as e: print('Error: {}'.format(e))
        self.stop_spin()

    # option button commands ----------------------------------------------------

    def show_overview(self):
        tab = self.current_tab()
        if tab == self.HelpTab:
            self.show_welcome()
        elif tab == self.NormBallTab:
            self.show_nb_overview()
        elif tab == self.DualNormBallTab:
            self.show_dnb_overview()

    def show_welcome(self):
        datadir = os.path.dirname(__file__)
        path = os.path.join(datadir,'data','welcome.txt')
        with open(path, 'r') as file:
            self.welcome = tk.Text(self.HelpTab, wrap='none', bg="#DADADA", fg="#444444", height=34, width=103, borderwidth=-2, padx=-1, pady=-1, insertwidth=0)
            self.welcome.grid(row=0, column=0, sticky='nswe')
            for line in file.readlines():
                self.welcome.insert(tk.END, line)

    def show_nb_overview(self):
        generate_nb_overview(self)

    def show_dnb_overview(self):
        generate_dnb_overview(self)



    def show_summary(self):
        computing=True
        t = threading.Thread(target=generate_summary,args=(self,))
        t.start()
        self.start_spin()

    def show_vertices(self):
        tab = self.current_tab()
        self.clear_tab(tab)
        computing = True
        if tab == self.NormBallTab:
            t = threading.Thread(target=make_all_vertices_table, args=(self,tab))
        else:
            t = threading.Thread(target=make_dnb_vertices_table, args=(self,tab))
        t.start()
        self.start_spin()


    def show_faces(self):
        P = self.ball.polyhedron()
        dim = P.dim()-1
        make_facets_table(self, self.NormBallTab, dim)


    def show_facets(self):
        computing = True
        self.clear_tab(self.NormBallTab)
        dim = int(self.DimEntry.get())
        t = threading.Thread(target=make_facets_table, args=(self, self.NormBallTab, dim))
        t.start()
        self.start_spin()


    def show_all_qtons(self):
        computing = True
        self.clear_tab(self.QtonsTab)
        t = threading.Thread(target=make_qtons_table,args=(self,))
        t.start()
        self.start_spin()

    def show_hasse(self):
        computing = True
        tab = self.current_tab()
        self.clear_tab(tab)
        t = threading.Thread(target=make_hasse, args=(self, tab))
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
        online = False
        tab = self.current_tab()
        self.clear_tab(tab)
        t = threading.Thread(target=_show_polyhedron, args=(self, online, tab))
        t.start()
        self.start_spin()


#    def cancel_button_click(self):
#        self.parent.destroy()


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
    icon_path = os.path.join(datadir, 'images', 'icon', 'icon.png')
    img = tk.PhotoImage(file=icon_path)
    root.tk.call('wm', 'iconphoto', root._w, img)
    root.mainloop()



if __name__ == '__main__':
    start()

