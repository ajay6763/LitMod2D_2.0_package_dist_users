#####################################################################################################
# This is the main python program of the gui 
#####################################################################################################
import Tkinter as tk
from Tkinter import *
import ttk
LARGE_FONT=("Times",12)
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry.polygon import  Polygon as Shape
from shapely.geometry import  Point, LinearRing
from matplotlib.widgets import Button, Cursor
from matplotlib.figure import Figure
from traits.api import HasTraits, Str, Int, Float, Enum, Array
from traitsui.api import View, Item, Group, HSplit, Handler 
from traitsui.menu import OKButton, CancelButton, ApplyButton,UndoButton

from time import sleep, strftime
from traits.api import *
from matplotlib.widgets import Button
import matplotlib.animation as animation
import math
import sys, string, os
import matplotlib.image as mpimg

import webbrowser # for web links
#### set plot figure before hand
f = Figure(figsize=(5,5),dpi=100)
a = f.add_subplot(111)

plt.ioff()
#from mpl_figure_editor import MPLFigureEditor
import my_lib
#from matplotlib.widgets import Button
#from matplotlib import style
#style.use("classic") #### style o

# getting the path of the Litmod program
litmod_path=os.getcwd()


class Model_setup(HasTraits):
	X_length = 1000
	Y_length = 400
	X_resolution = 10
	# This decides what to be shown in the GUI
	view1 = View(Group(Item(name = 'X_length'),
				 Item(name = 'Y_length'),
				 Item(name = 'X_resolution'),
				 show_border = True),
				 title="Model Setup",
				 buttons = [OKButton, CancelButton])

class Obs_data(HasTraits):
	topo = "test_top.dat"
	bouger = "test_bou.dat"
	geoid = "test_geo.dat"
	thermo_data = 5
	grain_size = 5
	oscillation = 75
	anomaly = 0
	print_option = 1
	# This decides what to be shown in the GUI
	view1 = View(Group(Item(name = 'topo'),
				 Item(name = 'bouger'),
				 Item(name = 'geoid'),
				 Item(name = 'thermo_data'),
				 Item(name = 'grain_size'),
				 Item(name = 'oscillation'),
				 Item(name = 'anomaly'),
				 Item(name = 'print_option'),
				 show_border = True),
				 title="Observed data",
				 buttons = [OKButton, CancelButton])




##################################################################################################################3
################ main program
##################################################################################################################3

class lithmod(tk.Tk):
	def __init__(self, master=None, *args, **kwargs):
		tk.Tk.__init__(self, master, *args, **kwargs) # initialise tkinter
		
		#tk.Tk.wm_iconbitmap(self,bitmap ="@gdl.ico") ### to add .wm_iconbitmap(bitmap = "myicon.ico")the logo in the
		#.wm_iconbitmap(bitmap = "logo.ico")
		tk.Tk.wm_title(self,"LitMod")
		container=tk.Frame(self)
		container.pack(side="top",fill="both",expand=True)		
		container.grid_rowconfigure(0,weight=1)
		container.grid_columnconfigure(0,weight=1)
    
		### specify dictionary
		self.frames={} ## all frames

		for F in (StartPage, Help):
			frame=F(container,self)
			self.frames[F]=frame
			frame.grid(row=0,column=0,sticky="nsew")
		
		self.show_frame(StartPage)  ## which frame to show_frame

	def show_frame(self,cont):
		frame=self.frames[cont]
		frame.tkraise() ## tkinter function to raise the frame to up
	def web_link(self):
		webbrowser.open_new("https://sites.google.com/site/ictjagdl/people")
    	#return webbrowser.open_new("http://www.google.com")


		
#### Add example StartPage
class StartPage(tk.Frame):
	
	def __init__(self, parent, controller):
		tk.Frame.__init__(self,parent)
		label=tk.Label(self,text="Welcome to LitMod",font=LARGE_FONT)
		label.pack(pady=10,padx=10)
		button2=ttk.Button(self,text="Build Model",
	 		command=lambda:my_lib.build_model(litmod_path))
		button2.pack()
		button4=ttk.Button(self,text="Load Model",
	 		command=lambda:my_lib.load_model(litmod_path))
		button4.pack()
		button1=ttk.Button(self,text="Help!",
	 		command=lambda:controller.show_frame(Help))
		button1.pack()
		button3=ttk.Button(self,text="Tips! \n \
	 	  1.Build \n * When entering the properties ofthe bodies please keep the format in which default values are.\
	 		\n*To draw bodies you should have a sketch of the profile where you should know where bodies connect. \n * If you want to know more in property editor fields, just hold the cursor in the field a balck box will appear with info.\n\
	 		2.Load Model\n* To load a model you should have a litmod.inp and bodies.out file. If you only have litmod.inp then just run LitMod it will produce bodies.out\n or some file with infor of node points just copy that to the bodies.out.\n \
	 		3. Input,output files \n * Put all input files in one folder where you should also have thermodynamic tables from Generator.\n *Ouput files name are already entered. You can use them to plot some where else(e.g. GMT) ")
		button3.pack()
	 	

class Help(tk.Frame):

	def __init__(self,parent,controller):
		tk.Frame.__init__(self,parent)
		label=tk.Label(self,text="Help",font=LARGE_FONT)
		label.pack(pady=10,padx=10)
		button2=ttk.Button(self,text="LitMod webpage",
	 		command=lambda:controller.web_link())
	 	#button2.bind(web_link(self.event))
		button2.pack()
		button1=ttk.Button(self,text="Back to home",
	 		command=lambda:controller.show_frame(StartPage))
		button1.pack()
	 	
	 	#button2=ttk.Button(self,text="Page2",
	 	#	command=lambda:controller.show_frame(Plot))
	 	#button2.pack()
"""
class Plot(tk.Frame):

	def __init__(self,parent,controller):
		tk.Frame.__init__(self,parent)
		label=ttk.Label(self,text="Cross-section",font=LARGE_FONT)
		label.pack(pady=10,padx=10)
                
		button1=ttk.Button(self,text="Load Body",
	 		command=lambda:load_body("body1x.txt","body1y.txt"))
	 	button1.pack()
	 	
	 	button2=ttk.Button(self,text="Back to home",
	 		command=lambda:controller.show_frame(StartPage))
	 	button2.pack()
	 	button3=ttk.Button(self,text="Help",
	 		command=lambda:controller.show_frame(Help))
	 	button3.pack()


class Page3(tk.Frame):

	def __init__(self,parent,controller):
		tk.Frame.__init__(self,parent)
		label=ttk.Label(self,text="Add body",font=LARGE_FONT)
		label.pack(pady=10,padx=10)

	 	button1=ttk.Button(self,text="Back to home",
	 		command=lambda:controller.show_frame(StartPage))
	 	button1.pack()

	 	## plot data it is in the background
	 	## now get the plot to the tkinter canvas
	 	canvas=FigureCanvasTkAgg(f,self)
	 	canvas.show()
	 	canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)

	 	## add the matplotlib toolbar
	 	toolbar = NavigationToolbar2TkAgg(canvas, self) # adding 
	 	toolbar.update()
	 	canvas._tkcanvas.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
"""

if __name__ == "__main__":
	app=lithmod()
	app.mainloop()
