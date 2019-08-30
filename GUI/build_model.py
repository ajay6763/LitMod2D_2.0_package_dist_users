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
from time import sleep
from traits.api import *
from scipy import *
from matplotlib.widgets import Button
import matplotlib.animation as animation
import math
import sys, string, os
import matplotlib.image as mpimg

plt.ioff()
#from mpl_figure_editor import MPLFigureEditor
from scipy import *
import my_lib
#from matplotlib.widgets import Button
from matplotlib import style
style.use("dark_background") #### style o

class Welcome(HasTraits):
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


model=Welcome()
model.configure_traits()


#crust = Prop()
body=[]
body.append(my_lib.Body(1))


my_lib.draw_line(model.X_length,model.Y_length,model.X_resolution)
"""
####################################################################33
#########3 writing litmod.inp file
f=open("litmod.inp","w")
obs=Obs_data()
obs.configure_traits(view="view1")
f.writelines('{:74s} {:s}\n'.format(obs.topo,"ELEVin"))
f.writelines('{:74s} {:s}\n'.format(obs.bouger,"Bouguer"))
f.writelines('{:74s} {:s}\n'.format("      ","Free-air"))
f.writelines('{:74s} {:s}\n'.format(obs.geoid,"GEOIDin"))
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
f.writelines(["%s                                                                        %s\n" % (obs.thermo_data,"Thermo D")])
f.writelines(["%s" % ("comment\n")])
#### writing some info about the model
f.writelines(["   %s    %s    %s    %s   %s   %s\n" % (len(np.unique(body[0].material)),len(body[0].x),obs.anomaly,obs.grain_size,obs.oscillation,obs.print_option)])
#### writing body property info
for i in range(len(body[0].x)):
	print body[0].x[i]
	print body[0].name[i]
	print body[0].color[i]
	f.writelines(["      %s      %s      %s\n" % (i+1,body[0].material[i],body[0].name[i],)])
	f.writelines(["  %s  %s  %s  %s  %s  %s  %s\n" % (body[0].heat_produc[i],body[0].heat_produc_depth[i], \
						body[0].conducX[i],body[0].conducY[i],"0.000E+00","0.000E+00","0.000E+00")])
	print "printing density depth"
	f.writelines(["  %s  %s  %s  %s\n" % (body[0].density[i],"0.000E+00","0","0.000E+00")])
###### writing body cordinates	
for j in range(len(body[0].x)):
	f.writelines(["%s                 %s\n" % (j+1,body[0].name[j])])
	for i in range(len(body[0].x[j])):
		if i==(len(body[0].x[j])-1) or (i+1)%4 ==0:
			f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(body[0].x[j][i]*1000)),".",int(float(body[0].y[j][i]*1000)),"."))
			f.writelines("\n")
		else: 
			f.writelines('{:7d}{:s}{:7d}{:s}'.format(int(float(body[0].x[j][i]*1000)),".",int(float(body[0].y[j][i]*1000)),"."))
#### writing the boundary conditions for temprature
f.writelines(["%s %s %s %s %s\n" % ("0","1","1","0","0")])
#### writing the boundary conditions for heat flow
f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","0.,","0.,","0.")])
#### writing the boundary conditions value for temerature
f.writelines(["%s%s%s%s%s\n" % ("0.,","0.,","1320.,","0.,","0.")])
### writing nodes information
f.writelines(["%s " " %s " " %s " " %s\n" % ("0.",float(str(model.X_resolution))*1000,len(range(0,int(str(model.X_length)),int(str(model.X_resolution)))),"96")])
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
#print body[0].name[1]
print len(body[0].x)
"""