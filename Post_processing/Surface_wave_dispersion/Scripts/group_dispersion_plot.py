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

fig  = plt.figure(figsize=(6,8 ))
dist=[100,200,250,300,400,600]
ax_group = fig.add_subplot(121)
#dist=[8,27,48,67,89,99]
for i in range(len(dist)):
	s=str(str(dist[i])+"km_R_U_SURF96.inp")
	data=np.loadtxt(s,usecols=(5,6,7))
	ax_group.plot(data[:,1],data[:,0],label=str(dist[i])+"km")
#ax_ray.grid(True)
ax_group.invert_yaxis()
ax_group.grid(True,linestyle='--',color='black',linewidth=0.15)
plt.xticks(fontsize=axes_label_font_size)
plt.yticks(fontsize=axes_label_font_size)
plt.ylabel('Time period ($s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.xlabel('Group velocity ($km/s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.legend( fancybox=True,shadow=True, framealpha=0.5,loc='lower right',fontsize=legend_font_size)

ax_phase = fig.add_subplot(122)
#dist=[8,27,48,67,89,99]
for i in range(len(dist)):
	s=str(str(dist[i])+"km_R_C_SURF96.inp")
	data=np.loadtxt(s,usecols=(5,6,7))
	ax_phase.plot(data[:,1],data[:,0],label=str(dist[i])+"km")
#ax_ray.grid(True)
ax_phase.invert_yaxis()
ax_phase.grid(True,linestyle='--',color='black',linewidth=0.15)
plt.yticks(fontsize=axes_label_font_size)
plt.xticks(fontsize=axes_label_font_size)
plt.ylabel('Time period ($s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.xlabel('Phase velocity ($km/s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.legend( fancybox=True,shadow=True, framealpha=0.5,loc='lower right',fontsize=legend_font_size)


plt.savefig('Rayleigh_Group_Phase_dispersion.jpg', dpi=400)  
plt.show()
