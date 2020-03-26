#####################################################################################################
# This is the main python program of the gui 
#####################################################################################################
import Tkinter as tk
import ttk
from matplotlib.figure import Figure
from traitsui.menu import OKButton, CancelButton, ApplyButton,UndoButton
import os
import webbrowser # for web links
#############################
### This the liberary where all the GUI functions resides
import my_lib

############################
### getting the path of the Litmod program
litmod_path=os.getcwd()
LARGE_FONT=("Times",12)

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
################################### 
### Start page of the GUI
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
################################### 
### Help page of the GUI
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
if __name__ == "__main__":
	app=lithmod()
	app.mainloop()
