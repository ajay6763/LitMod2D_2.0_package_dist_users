#### ploting liberary
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import os
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
dist=[30,60,80,95]

for i in range(len(dist)):

	s=str(str(dist[i])+"km.xyz")
	data=np.loadtxt(s)
	#data[:,1]=data[:,1]/max(data[:,1])
	ax_ray.stackplot(data[101::,1]+i,data[101::,0],labels=str(dist[i]*1000))
ax_ray.grid(True)
#ax_ray.invert_yaxis()
#plt.gca().invert_yaxis()
ax_ray.grid(True,linestyle='--',color='black',linewidth=0.15)

#plt.ylim(-40,0)
plt.xticks(fontsize=axes_label_font_size)
plt.yticks(fontsize=axes_label_font_size)
plt.xlabel('RF apmlitude',fontsize=axes_label_font_size, fontweight='bold')
plt.ylabel('Time ($s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.legend( fancybox=True,shadow=True, framealpha=0.5,loc='lower right',fontsize=legend_font_size)
plt.savefig('RF.jpg', dpi=400)  
plt.show()
