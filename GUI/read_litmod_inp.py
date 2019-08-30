import wx
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from shapely.geometry.polygon import  Polygon as Shape
#from shapely.geometry import Polygon, Point, LinearRing
from shapely.geometry import  Point, LinearRing
from matplotlib.widgets import Button
from matplotlib.figure import Figure
from traits.api import HasTraits, Str, Int, Float, Enum
from traitsui.api import View, Item, Group
from traitsui.wx.themed_text_editor import \
    ThemedTextEditor # for different themes in the menu
from traitsui.menu import OKButton, CancelButton, ApplyButton,UndoButton
import traitsui
from threading import Thread
from time import sleep
from traits.api import *
from traitsui.api import View, Item, Group, HSplit, Handler
from traitsui.menu import NoButtons
#from mpl_figure_editor import MPLFigureEditor
from scipy import *
from matplotlib.widgets import Button
import matplotlib.animation as animation
import math
from matplotlib.widgets import Cursor # for fancy cursor
import plot_lib

import sys, string, os

from matplotlib import style
style.use("dark_background") #### style o

class Canvasread(HasTraits,object):
    name = "mantle"
    color = Enum("red", "blue", "green","yellow","orange","pink","brown","purple",)
    ref_density = "3.245E+03"
    Body_number= "0"
    Depth_density_factor= "0.000E+00"
    Temp_density_factor= "0"
    Thermal_expansion= "0.000E+00"

    material = 99
    ### thermal conductivity parameters
    Kx ="2.100E+00"
    Ky ="2.100E+00"
    Kz ="0.000E+00" #angle of anisotropy
    Temp_conduc_factor="0.000E+00"
    Depth_conduc_factor="0.000E+00"
    ### heat production parameters
    radiogenic_heat_producion = "2.000E-08"
    depth_heat_factor="0.000E+00"

    #shape=Shape([(0,0),(1000,0),(1000,-400),(0,-400)])
    shape=Shape([(0,0),(1000,0),(1000,-400),(0,-400)])
    x_s = "0.0"
    y_s = "-100"
    x_e = "1000"
    y_e = "-100"
    body_type= 1
    #############################################################################3
    #### selferved data and some parameteres for input file
    topo = File #"test_top.dat"
    bouger = File #"test_bou.dat"
    geoid = File #"test_geo.dat"
    thermo_data = 5
    grain_size = 5
    oscillation = 75
    anomaly = 0
    print_option = 1
    topo_response= Enum("no","yes")
    # This decides what to be shown in the GUI
    
    view1 = View(Group(Item(name = 'name'),
        Item(name = 'Body_number'),
        Item(name = 'color', label = "click to select"),
             Item(name = 'material'),
             Group(Item(name = 'ref_density',label = "Reference Density(kg/m3)"),Item(name = 'Depth_density_factor',label="P dependence factor"),\
                Item(name = 'Temp_density_factor',label="T dependence factor"),Item(name = 'Thermal_expansion',label="Thermal expansion  factor"),\
                orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Density Parameters"),
             Group(Item(name = 'Kx',label="Conductivity in X"),Item(name = 'Ky',label="Conductivity in Y"),\
                Item(name = 'Temp_conduc_factor',label="T Conductivity Factor"),Item(name = 'Depth_conduc_factor',label="T Conductivity Factor"),\
                Item(name = 'Kz',label="Anisotropy angle"),orientation='horizontal', layout='tabbed', springy=True,show_labels=True,label="Thermal Conductivity"),
             Group(Item(name = 'radiogenic_heat_producion',label="Heat production A0(W\m3)"),Item(name = 'depth_heat_factor',label="Reference thickness of exponential production(m)"),\
                orientation='horizontal', layout='tabbed', springy=True,show_labels=True, label="Heat Production  Parameters"),
             #Group(Item(name = 'radiogenic_heat_producion')),
             Item(name = 'body_type'),
             Item(name = 'x_s'),
             Item(name = 'y_s'),
             Item(name = 'x_e'),
             Item(name = 'y_e'),
             show_border = True),
             title="Property Editor",
             buttons = [OKButton, CancelButton]) 
    view2 = View(Group(Item(name = 'topo'),
                 Item(name = 'bouger'),
                 Item(name = 'geoid'),
                 Item(name = 'thermo_data'),
                 Item(name = 'grain_size'),
                 Item(name = 'oscillation'),
                 Item(name = 'anomaly'),
                 Item(name = 'print_option'),
                 show_border = True),
                 title="selferved data",
                 buttons = [OKButton, CancelButton])   
    view3 = View(Group(Item(name = 'Body_number',label="Enter number of the Body to edit"),
             show_border = True),
             title="Property Editor",
             buttons = [OKButton, CancelButton])


    def __init__(self,ax,save_option,add_option,edit_option,plot_option,run_option,X_length,Y_length,X_resolution,\
    	l,b,material,name,density,density_depth,density_temp,thermal_expan_factor,conducX,conducY,conducZ,\
		Depth_conduc_factor,Temp_conduc_factor,heat_produc,heat_produc_depth,n_bod):
        self.save_option = save_option
        self.add_option = add_option
        self.edit_option = edit_option
        self.plot_option = plot_option
        self.run_option = run_option
        self.X=X_length
        self.Y=float(str(Y_length))
        self.reso=float(str(X_resolution))
        self.Z=[0,  -0.5,  -1.5,  -2,  -2.5,  -3,  -3.5,  -4, \
        -4.5,  -5.5,  -6,  -6.5,  -7,  -8,  -9, -10, \
        -11, -12, -13, -14, -15, -17, -19, -21, \
        -23, -25, -27, -29, -31, -33, -35, -37, \
        -40, -43, -46, -49, -52, -55, -58, -61, \
        -64, -67, -70, -75, -80, -85, -90, -95, \
        -100,-105,-110,-115,-120,-125,-130,-135, \
        -140,-145,-150,-155,-160,-165,-170,-175, \
        -180,-185,-190,-195,-200,-205,-210,-215, \
        -220,-225,-230,-240,-250,-260,-270,-280, \
        -290,-300,-305,-310,-315,-320,-330,-340, \
        -350,-360,-370,-380,-390,-400,-410,-420]
        self.ax = ax
        self.ax1 =ax
        self.n = name
        self.n_bod=n_bod
        self.c = []
        self.mat = material
        temp_color=["green","blue",'yellow','red','pink','orange','purple','white','black']
        for i in range(int(str(self.n_bod))):
        	self.c.append(temp_color[i])

        # density parameters
        self.dens = density
        self.dens_depth =density_depth
        self.dens_temp = density_temp
        self.therm_expan = thermal_expan_factor
        # conductivity parameters
        self.conducx=conducX
        self.conducy=conducY
        self.conducz=conducZ
        self.depth_conduc_factor=Depth_conduc_factor
        self.temp_conduc_factor=Temp_conduc_factor

        self.heat_produ = heat_produc
        self.heat_produ_depth = heat_produc_depth
        
        
        self.b = b
        self.x_st= []
        self.y_st= []
        self.x = []
        self.y = []
        self.x_end = [] 
        self.y_end = []
        self.l = l
        # Set limits to unit square
        #self.ax.set_xlim((0,1000))
        #self.ax.set_ylim((-400,0))
        #ax.set_title('Profile')
        #### to add image in the background
        #img = plt.imread("figure_2.png")
        #self.ax.imshow(img)

        self.ax.set_xlim((0, self.X))
        self.ax.set_ylim((-self.Y,0))
        major_xticks = np.arange(0, self.X, 50)
        major_yticks = np.arange(-self.Y, 0, 25)                                              
        self.ax.set_yticks(major_yticks)
        self.ax.set_xticks(major_xticks)  
        self.ax.grid(True)
        self.ax.set_title('LEFT: new point, MIDDLE: delete last point, RIGHT: close polygon')
        #self.mouse_button = {1: self._add_point, 2: self._delete_point, 3: self._close_polygon}
    def plot_model(self,event):
     	print "herre tou dumvc"
        if event.inaxes == self.plot_option:
            #temp_color=['green','blue','yellow','red','pink','orange','purple','white','black']
            for i in range(len(self.l)):
            	print "went inside"
            	print self.l[i]
                poly = Polygon(list(zip(self.l[i],self.b[i])),facecolor=self.c[i],label=self.mat[i], lw=1)
                
                #self.ax.text(self.l[i][len(self.l[i])/2], self.b[i][len(self.b[i])/2], "asdf", fontsize=15,color="black")
                self.ax.add_patch(poly)
                self.l_p=(max(self.l[i][:]) + min(self.l[i][:]))/2
                self.b_p=(max(self.b[i][:]) + min(self.b[i][:]))/2 
                self.ax.text(self.l_p, self.b_p, i+1, fontsize=15,color="black")
                #plt.draw()
            plt.draw()    
    def run_model(self,event):
    	print "running the model"
        if event.inaxes == self.run_option:
            print "running the model"
            os.system("./Litmod4.0")
            plot_lib.plot_lib(len(range(0,int(float(self.X)),int(float(self.reso))))+1,self.l,self.b)
    def edit_body(self,event):
        # If the mouse pointer is not on the canvas, ignore buttons
        if event.inaxes == self.edit_option:
            print "in edit"
            #n_b=input("Enter index of the body to edit ")
            self.configure_traits(view='view3')
            n_b=int(self.Body_number)
            self.name = self.n[n_b-1]
            self.color = self.c[n_b-1]

            self.ref_density = self.dens[n_b-1]
            self.Depth_density_factor= self.dens_depth[n_b-1]
            self.Temp_density_factor=self.dens_temp[n_b-1]
            self.Thermal_expansion=self.therm_expan[n_b-1]
            self.material = self.mat[n_b-1]
            
            self.Kx =self.conducx[n_b-1]
            self.Ky =self.conducy[n_b-1]
            self.Kz =self.conducz[n_b-1]
            self.Depth_conduc_factor=self.depth_conduc_factor[n_b-1]
            self.Temp_conduc_factor=self.temp_conduc_factor[n_b-1]

            self.radiogenic_heat_producion =self.heat_produ[n_b-1]
            self.depth_heat_factor=self.heat_produ_depth[n_b-1]
            """
            self.x_s = self.x_st[n_b-1]
            self.y_s = self.y_st[n_b-1]
            self.x_e = self.x_end[n_b-1]
            self.y_e = self.y_end
            self.body_type= 1
            """
            self.configure_traits(view='view1')
            #n_b=int(str(self.Body_number))
            #x_s=input("Enter X start point of body: ")
            #y_s=input("Enter Y start point of body: ")
            #print n_bod-1
            print self.l[n_b-1]
            print self.b[n_b-1]
            self.n[n_b-1]=self.name
            self.c[n_b-1]=self.color
            self.conducx[n_b-1]=self.Kx
            self.conducy[n_b-1]=self.Ky
            self.conducz[n_b-1]=self.Kz
            self.depth_conduc_factor[n_b-1]=self.Depth_conduc_factor
            self.temp_conduc_factor[n_b-1]=self.Temp_conduc_factor
            self.mat[n_b-1]=self.material
            self.heat_produ[n_b-1]=self.radiogenic_heat_producion
            self.heat_produ_depth[n_b-1]=self.depth_heat_factor
            self.dens[n_b-1]=self.ref_density
            self.dens_depth[n_b-1]=self.Depth_density_factor
            self.dens_temp[n_b-1]=self.Temp_density_factor
            self.therm_expan[n_b-1]=self.Thermal_expansion
        
            #for i in range(len(self.l)):
            #    print i
            #self.path.set_data(self.l[n_b-1],self.b[n_b-1])
            self.ax.plot(self.l[n_b-1],self.b[n_b-1],color=self.c[n_b-1],lw=1,marker='o')
            poly = Polygon(list(zip(self.l[n_b-1],self.b[n_b-1])),facecolor=self.c[n_b-1],label=self.mat[n_b-1], lw=1)
            self.ax.add_patch(poly)
            poly=[]
            #poly = Polygon(list(zip(self.l[self.n_bod-1], self.b[self.n_bod-1])),facecolor=self.color, lw=1,animated=True)
            #self.ax.update.add_patch(poly)
            #self.ax.add_patch(poly)
            plt.draw()
    def save_body(self,event):
        # If the mouse pointer is not on the canvas, ignore buttons
        if event.inaxes == self.save_option:
            #Prop.configure_tself.density[i]raits(View=self.view1)
            self.configure_traits(view='view2')

            #### writing bodies coordinates ina file which is to be read if you want to read this inpout file again
            f=open("bodies.out","w")
            f.writelines(">\n")
            for j in range(len(self.l)):
                for i in range(len(self.l[j])):
                    if i==(len(self.l[j])-len(self.l[j])):
                        f.writelines(">\n")
                    f.writelines(["%f \t %f \n" % (self.l[j][i] , self.b[j][i])])
            f.close()
            ####################################################################33
            #########3 writing litmod.inp file
            f=open("litmod.inp","w")
            f.writelines('{:74s} {:s}\n'.format(self.topo,"ELEVin"))
            f.writelines('{:74s} {:s}\n'.format(self.bouger,"Bouguer"))
            f.writelines('{:74s} {:s}\n'.format("      ","Free-air"))
            f.writelines('{:74s} {:s}\n'.format(self.geoid,"GEOIDin"))
            f.writelines('{:74s} {:s}\n'.format("topoout.dat","Topograp"))
            f.writelines('{:74s} {:s}\n'.format("bouguer_out.data","Bouguer"))
            f.writelines('{:74s} {:s}\n'.format("           ","Free air"))
            f.writelines('{:74s} {:s}\n'.format("        ","Mixed"))
            f.writelines('{:74s} {:s}\n'.format("geoid2D.dat","GRAVout"))
            f.writelines('{:74s} {:s}\n'.format("tempout.dat","Temperat"))
            f.writelines('{:74s} {:s}\n'.format("fluxout.dat","SHF out"))
            f.writelines('{:74s} {:s}\n'.format("           ","CSTR out"))
            f.writelines('{:74s} {:s}\n'.format("           ","ESTR out"))
            f.writelines('{:74s} {:s}\n'.format("           ","COL out"))
            f.writelines('{:74s} {:s}\n'.format("           ","BOUN out"))
            f.writelines('{:74s} {:s}\n'.format("bodies.dat","BODY out"))
            f.writelines('{:74s} {:s}\n'.format("           ","ELEM out"))
            f.writelines('{:74s} {:s}\n'.format("           ","ISTR out"))
            f.writelines('{:74s} {:s}\n'.format("           ","P/T out"))
            f.writelines('{:74s} {:s}\n'.format("dens.dat","Density"))
            f.writelines('{:74s} {:s}\n'.format("           ","RECA out"))
            f.writelines(["%s                                                                        %s\n" % (self.thermo_data,"Thermo D")])
            f.writelines(["%s" % ("comment\n")])
            #### writing some info about the model
            f.writelines(["   %s    %s    %s    %s   %s   %s\n" % (len((self.mat)),len(self.l),self.anomaly,self.grain_size,self.oscillation,self.print_option)])
            #### writing body property info
            for i in range(len(self.l)):
                print self.l[i]
                print self.n[i]
                print self.c[i]
                f.writelines(["      %s      %s      %s\n" % (i+1,self.mat[i],self.n[i],)])
                f.writelines(["  %s  %s  %s  %s  %s  %s  %s\n" % (self.heat_produ[i],self.heat_produ_depth[i], \
                        self.conducx[i],self.conducy[i],self.conducz[i],self.temp_conduc_factor[i],self.depth_conduc_factor[i])])
                print "printing density depth"
                f.writelines(["  %s  %s  %s  %s\n" % (self.dens[i],self.dens_depth[i],self.dens_temp[i],self.therm_expan[i])])
            ###### writing body cordinates  
            for j in range(len(self.l)):
                f.writelines(["%s                 %s\n" % (j+1,self.n[j])])
                for i in range(len(self.l[j])):
                    if i==(len(self.l[j])-1) or (i+1)%4 ==0:
                        f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
                        f.writelines("\n")
                    else: 
                        f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(self.l[j][i]*1000)),".",int(float(self.b[j][i]*1000)),"."))
            #### writing the boundary conditions for temprature
            f.writelines(["%s %s %s %s %s\n" % ("0","1","1","0","0")])
            #### writing the boundary conditions for heat flow
            f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","0.,","0.,","0.")])
            #### writing the boundary conditions value for temerature
            f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","1320.,","0.,","0.")])
            ### writing nodes information
            f.writelines(["%s " " %s " " %s " " %s\n" % ("0.",self.reso*1000,len(range(0,int(float(self.X)),int(float(self.reso))))+1,"96")])
            f.write("      0.   -500.  -1500.  -2000.  -2500.  -3000.  -3500.  -4000.\n\
  -4500.  -5500.  -6000.  -6500.  -7000.  -8000.  -9000. -10000.\n\
 -11000. -12000. -13000. -14000. -15000. -17000. -19000. -21000.\n\
 -23000. -25000. -27000. -29000. -31000. -33000. -35000. -37000.\n\
 -40000. -43000. -46000. -49000. -52000. -55000. -58000. -61000.\n\
 -64000. -67000. -70000. -75000. -80000. -85000. -90000. -95000.\n\
-100000.-105000.-110000.-115000.-120000.-125000.-130000.-135000.\n\
-140000.-145000.-150000.-155000.-160000.-165000.-170000.-175000.\n\
-180000.-185000.-190000.-195000.-200000.-205000.-210000.-215000.\n\
-220000.-225000.-230000.-240000.-250000.-260000.-270000.-280000.\n\
-290000.-300000.-305000.-310000.-315000.-320000.-330000.-340000.\n\
-350000.-360000.-370000.-380000.-390000.-400000.-410000.-420000.\n\
0,4,5,3,4\n\
0,0,0,0,0")   

##########################################################################################
############### draw_line use canvas to drwa body
##########################################################################################



class Body:
	def __init__(self,number):
		self.number = []
		self.name =[]
		self.color =[]
		self.x = []
		self.y = []
		self.material = []
		# density parameters
		self.density = []
		self.density_depth = []
		self.density_temp = []
		self.thermal_expan = []
		# conductivity parameters
		self.conducX=[]
		self.conducY=[]
		self.conducZ=[]
		self.Depth_conduc_factor=[]
		self.Temp_conduc_factor=[]
		### heat production parameters
		self.heat_produc = []
		self.heat_produc_depth = []
##################################################################################################
# reading general info 
##################################################################################################

fopen = open("litmod.inp")
content=fopen.readlines()
fopen.close()

## reading general info about the model setup
n_mat=content[23].split()[0]
n_bod=content[23].split()[1]
iel_c=content[23].split()[2]
dsize=content[23].split()[3]
iospe=content[23].split()[4]
iprint=content[23].split()[5]

# making array of body ojbjects
#bodies = []

#for i in range(int(str(n_bod))+1):
#	bodies.append(Body(i))
body=[]
body.append(Body(1))
##################################################################################################
# reading bodies properties
##################################################################################################

j=1
for i in range(int(str(n_bod))):
    # here appending the body
    body[0].name.append([])
    body[0].number.append([])
    body[0].material.append([])
    #body[0].color.append([])
    #density parameters
    body[0].density.append([])
    body[0].density_depth.append([])
    body[0].density_temp.append([])
    body[0].thermal_expan.append([])
    # conductivity parameters    
    body[0].conducX.append([])
    body[0].conducY.append([])
    body[0].conducZ.append([])
    body[0].Depth_conduc_factor.append([])
    body[0].Temp_conduc_factor.append([])
    # heat production parameter
    body[0].heat_produc.append([])
    body[0].heat_produc_depth.append([])
    #here will store the info
    #for j in range(24+(3*i),24+(3*i),1)
    temp=content[24+(3*i)].split()
    print temp
    body[0].number[i]=int(str(temp[0]))
    body[0].material[i]=int(str(temp[1]))
    body[0].name[i]=(str(temp[2]))

    temp=content[24+(3*i)+1].split()
    # heat production parameters
    body[0].heat_produc[i]=(str(temp[0]))
    body[0].heat_produc_depth[i]=(str(temp[1]))
    # conductivity parameters
    body[0].conducX[i]=(str(temp[2]))
    body[0].conducY[i]=(str(temp[3]))
    body[0].conducZ[i]=(str(temp[4]))
    body[0].Temp_conduc_factor[i]=(str(temp[5]))
    body[0].Depth_conduc_factor[i]=(str(temp[6]))

    # density factors
    temp=content[24+(3*i)+2].split()
    body[0].density[i]=(str(temp[0]))
    body[0].density_depth[i]=(str(temp[1]))
    body[0].density_temp[i]=(str(temp[2]))
    body[0].thermal_expan[i]=(str(temp[3]))	

#for i in range(1,int(str(n_bod))+1):
print body[0].number
print body[0].material
print body[0].name
print body[0].heat_produc
print body[0].heat_produc_depth
print body[0].conducX
print body[0].conducY
print body[0].conducZ
print body[0].density
print body[0].density_depth
# using bodies.out
f=open("bodies.out")
coords=f.readlines()
j=-1
for i in range(1,len(coords)):
	
	temp=coords[i].split()
	#print temp
	if(temp[0]!='>'):
		#print i
		#print j
		body[0].x[j].append(float(str(temp[0])))
		body[0].y[j].append(float(str(temp[1])))
	if(temp[0]=='>'):
		j=j+1
		body[0].x.append([])
		body[0].y.append([])
		#print i
		#print j
print "printing bodies loaded"    
print body[0].x
print body[0].y
## reading info about model like length resolution and stuff
info=content[-15].split()
X_resolution=float(str(info[1]))/1000.0 #10
X_nodes=int(str(info[2]))
X_length=(X_resolution*X_nodes)-X_resolution  #1000
Y_length=400 



##################
fig = plt.figure(0)
fig.canvas.set_window_title('Litmod Read Model')
ax = fig.add_subplot(111)

#fig, ax = plt.subplots()
cursor = Cursor(ax, useblit=True, color='white', linewidth=0.5)
plt.subplots_adjust(bottom=0.2)
plt.xlabel('Distance (Km)')
plt.ylabel('Depth (Km)')
plot_option = plt.axes([0.4, 0.05, 0.1, 0.075])
run_option = plt.axes([0.51, 0.05, 0.1, 0.075])
save_option = plt.axes([0.62, 0.05, 0.1, 0.075])
add_option = plt.axes([0.73, 0.05, 0.1, 0.075])
edit_option = plt.axes([0.84, 0.05, 0.1, 0.075])
ax.set_title('Click and drag a point to move it')
ax.set_xlim((-2, 2))
ax.set_ylim((-2, 2))
cnv = Canvasread(ax,save_option,add_option,edit_option,plot_option,run_option,X_length,Y_length, \
	X_resolution,body[0].x,body[0].y,body[0].material,body[0].name, \
	body[0].density,body[0].density_depth,body[0].density_temp,\
	body[0].thermal_expan,body[0].conducX,body[0].conducY,body[0].conducZ,body[0].Temp_conduc_factor,body[0].Depth_conduc_factor,\
	body[0].heat_produc,body[0].heat_produc_depth,n_bod)
#print cnv.name
#print cnv.color
bnext = Button(add_option, 'add',color="black")
#bnext.on_clicked(cnv.add_body)
bprev = Button(save_option, 'save',color="black")
#bprev.on_clicked(crust.save)
bprev.on_clicked(cnv.save_body)
edit = Button(edit_option, 'edit',color="black")
edit.on_clicked(cnv.edit_body)
plot = Button(plot_option,'plot',color="black")
plot.on_clicked(cnv.plot_model)
run= Button(run_option,'run',color="black")
run.on_clicked(cnv.run_model)
        #bprev.on_clicked(crust.save)
#plt.connect('button_press_event',cnv.update_path)
#plt.connect('motion_notify_event',cnv.set_location)
plt.show()
        #