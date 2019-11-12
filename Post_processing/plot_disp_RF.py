#### ploting liberary
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib import gridspec
import numpy as np
import matplotlib as mpl
import os
from matplotlib.widgets import MultiCursor
from matplotlib import style
from pylab import *



##########################################
### Setting the paths to the data 
import os
cwd = os.getcwd()
disp_path=cwd+'/Surface_Wave_Dispersion_Curves/'
RF_path=cwd+'/Receiver_Functions/'
vel_path=cwd+'/Receiver_Functions/'
#import seaborn as sns
origin='lower'
label_size = 12
marker_size =2
fontsize = 12
legend_font_size=12
axes_label_font_size=12
linewidth=2
contour_font=10
#l=input('enter length:')
l=1000
#reso=input('enter resolution:')
reso=10

surf_type='R'
surf_type=input('Enter the Surface wave type. (Rayleigh-- "R", Love --"L" ):')
disp_type='U'
disp_type=input('Enter the Disperiosn type. (Phase-- "C" , Group-- "U" ): ')
max_depth = 400
min_depth = -10
x_nodes=l/reso
b=max_depth
hori=l/10
vert= (hori/l)*b
vert= (hori/l)*b

dist=[0,165,195,235,400,500,600,250,300,400,600]
dist = range(0,625,50)
dist = []
fig_len=len(dist)*2

fig, axes = plt.subplots(figsize=(fig_len,10),nrows=3, ncols=len(dist), sharex=False, sharey=False)
##########################################################
### Plotting Velocity
## Here deteminig the depth from CPS type velocity file
s=s=vel_path+str(str(dist[0])+"km_vel.dat") 
temp=np.loadtxt(s,skiprows=12)
depth=[]
for i in range(len(temp)):
  depth.append(np.sum(temp[0:i,0]))
#### now plottong the velocity profile
for i,d in enumerate(dist): #range(len(dist)):
  s=vel_path+str(str(d)+"km_vel.dat")
  data=np.loadtxt(s,skiprows=12)
  axes[0,i].plot(data[0:len(depth),2],depth[:],'k') #,label=str(dist[i])+"km")
  axes[0,i].grid(True,linestyle='--',color='black',linewidth=0.15)
  axes[0,i].invert_yaxis()
  axes[0,i].set_ylim([390,0])
  axes[0,i].set_xlim([2.5,4.9])
  axes[0,i].set_title(str(dist[i])+' km')
  axes[0,i].set_ylabel('Depth $km$)')
  axes[0,i].set_xlabel('Vs $km/s$)')
  
  #axes[0,i].set_xlim([5.0,4.8])

##########################################################
### Plotting RF
for  i,d in enumerate(dist):
  s=RF_path+str(str(d)+"km_vel.dat.2.5.eqr.xyz")
  data=np.loadtxt(s)
  axes[1,i].stackplot(data[:,1],data[:,0],color='k')
  try:
      s=RF_path+str(str(d)+"km.r.xyz")
      data=np.loadtxt(s)
      axes[1,i].plot(data[:,1],data[:,0],'r')  
  except:
    pass
  #axes[1,i].plot(data[:,1],data[:,0]) #,label=str(dist[i])+"km")
  axes[1,i].grid(True,linestyle='--',color='black',linewidth=0.15)
  axes[1,i].set_ylim([-1,20])
  axes[1,i].invert_yaxis()
  axes[1,i].set_ylabel('Time ($s$)')


##########################################################
### Plotting Disperiosn curves
for  i,d in enumerate(dist):
  try:
      s=disp_path+str(str(d)+"km_"+surf_type+"_"+disp_type+"_SURF96.inp")
      data=np.loadtxt(s,usecols=(5,6,7))
      axes[2,i].plot(data[:,1],data[:,0],'r')
  except:
    pass
  s=disp_path+str(str(d)+"km_"+surf_type+"_"+disp_type+"_SURF96.out")
  data=np.loadtxt(s,usecols=(5,6,7))
  axes[2,i].plot(data[:,1],data[:,0],'k')
  axes[2,i].invert_yaxis()
  axes[2,i].grid(True,linestyle='--',color='black',linewidth=0.15)
  axes[2,i].set_ylim((250,25))
  axes[2,i].set_xlim([3.0,5.0])
  axes[2,i].set_ylabel('Time period ($s$)')
  if disp_type=='U':
      axes[2,i].set_xlabel('Group velocity ($km/s$)')
  else:
      axes[2,i].set_xlabel('Phase velocity ($km/s$)')
  if surf_type=='R':
      axes[2,i].set_title('Rayleigh')
  else:
      axes[2,i].set_title('Love')
'''
#ax_ray.grid(True)
ax_group.invert_yaxis()
ax_group.grid(True,linestyle='--',color='black',linewidth=0.15)
'''


#plt.ylabel('Time period ($s$)',fontsize=axes_label_font_size, fontweight='bold')
#plt.xlabel('Group velocity ($km/s$)',fontsize=axes_label_font_size, fontweight='bold')
plt.tight_layout()
savefig('Vel_RF_DISP.jpg', dpi=400)
plt.show()