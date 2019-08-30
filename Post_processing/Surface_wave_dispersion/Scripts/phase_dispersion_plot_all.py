#### ploting liberary
import matplotlib as plt
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import numpy as np
import matplotlib as mpl
from traits.api import HasTraits, Str, Int, Float, Enum, Array, File, Directory
from traitsui.api import View, Item, Group, HSplit, Handler 
from traitsui.menu import OKButton, CancelButton, ApplyButton,UndoButton
import traitsui
import os
from matplotlib.widgets import MultiCursor
from matplotlib import style
import matplotlib.pyplot as plt
from matplotlib import gridspec

origin='lower'
label_size = 12
marker_size =2
fontsize = 12
legend_font_size=12
axes_label_font_size=12
linewidth=2
contour_font=10

fig  = plt.figure(figsize=(6,5 ))
ax_ray = fig.add_subplot(111)
dist=[20,60,85,95]
dist=[8,27,48,67,89,99]
dist=[10,25,30,35,40 ]# 15,20,25,30,35,40,45,50,55,60]
for i in range(len(dist)):
	s=str(str(dist[i])+"_SURF96.inp")
	data=np.loadtxt(s,usecols=(5,6,7))
	ax_ray.plot(data[:,0],data[:,1],label=str(dist[i]*10)+"km")
ax_ray.grid(True,linestyle='--',color='black',linewidth=0.15)

plt.xticks(fontsize=axes_label_font_size)
plt.yticks(fontsize=axes_label_font_size)
plt.xlabel('Time period ($s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.ylabel('Phase velocity ($km/s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.legend( fancybox=True,shadow=True, framealpha=0.5,loc='lower right',fontsize=legend_font_size)
plt.savefig('Rayleigh_Phase_dispersion.jpg', dpi=400)  
plt.show()
