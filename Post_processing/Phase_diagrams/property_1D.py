#### ploting liberary
import matplotlib as plt
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import numpy as np
#import matplotlib as mpl
import os
from matplotlib.widgets import MultiCursor
from matplotlib import style
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import *
import math
#import seaborn as sns
import scipy

cwd = os.getcwd()

origin='lower'
label_size = 10
marker_size =2
fontsize = 12
legend_font_size=10
axes_label_font_size=12
linewidth=2
contour_font=10

#########################################################
#### Anelastic attenuation model from Jackson et al., 2010
########################################################
AA1=816
alfa3=0.36
energi=293.0E03
volexp=1.20E-05
rgas=8.314472
pi=3.1415926
Qp=[]
Qs=[]
dsize=10
oscill=75

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='valid')
    return y_smooth

#######################
## functions to find the index
def lookup(T_LitMod,P_LitMod,T,P):
	index=[]
	dist=[]
#	dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
	dist=np.array((T[:]-T_LitMod)**2+(P[:]-P_LitMod)**2);
#	dist=np.array(((T_LitMod-T)**2+(P_LitMod-P)**2));
	index=dist.argmin(); 
	#print index, T_LitMod,P_LitMod
	return index
def find_mineral(mineral,mat_list,d):
	Ol=[]
	X=[]
	Z=[]
	d=d+1
	T_out=[]
	for j in range(len(mat_list)):
		s=str(str(mat_list[j])+"_profile")
		data=np.loadtxt(s)
		#T P Vp Vs Vp/Vs density wt% vol%
		s=str(str( mineral +"_"+ str(mat_list[j]))+"_FULL.txt")
		#print s
		try: 
			Ol_data=np.loadtxt(s)
			if len(Ol_data) ==0:
				print ('Length is Zero') 
				print s
				for i in range(len(data)):
					index=lookup(data[i,2],data[i,3],Ol_data[:,0],Ol_data[:,1])
					X.append(data[i,0])
					Z.append(-data[i,1])
					T_out.append(Ol_data[index,0])
					Ol.append(0)
				#pass
			else:
				print ('Length is not Zero') 
				print s
			
				## load the full tables for each material from Generator
				for i in range(len(data)):
					index=lookup(data[i,2],data[i,3],Ol_data[:,0],Ol_data[:,1])
					X.append(data[i,0])
					Z.append(-data[i,1])
					T_out.append(Ol_data[index,0])
					if d==2:
						parexp=math.exp((-(energi+(volexp*data[i,1])))/(rgas*(data[i,0]+273.0E0)))
						sqatt50=AA1*(((oscill*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3	
						cotp50= ((1.0E0/math.tan((pi*alfa3)/2.0E0))*sqatt50)*(2.0E0/9.0E0)
						Ol.append(Ol_data[index,d]*(1.0E0-cotp50))
					elif d==3:
						parexp=math.exp((-(energi+(volexp*data[i,1])))/(rgas*(data[i,0]+273.0E0)))
						sqatt50=AA1*(((oscill*(1.0E0/(dsize*1000.0E0)))*parexp))**alfa3	
						cots50= ((1.0E0/math.tan((pi*alfa3)/2.0E0))*sqatt50)*0.5E0
						Ol.append(Ol_data[index,d]*(1.0E0-cots50))
					else:
						Ol.append(Ol_data[index,d])
		except:
			pass	
					
	X=np.asarray(X)
	Z=np.asarray(Z)
	Ol=np.asarray(Ol)
	T_out=np.asarray(T_out)
	#### Make mesh grid
	#x,z=np.meshgrid(X, Z)
	#ol =np.tile(Ol,1)
	d=np.vstack((X.T,Z.T,Ol.T,T_out.T))
	#return X,Z,Ol,d.T
	#return x,z,ol
	return d.T
##############################
##### Main
#prop=input('Enter what property you want to plot: \n 1- P-wave velocities \n 2- S-wave velocities \n 3- Vp/Vs \n 4- Density \n 5- Weight % \n 6- Volume % \n 7- Mol% % \n ')

prop=input('Enter what property you want to plot: \n 1- P-wave velocities \n 2- S-wave velocities \n 3- Vp/Vs \n 4- Density \n 5- Weight % \n 6- Volume % \n  ')
prop_lw=input('Enter lower limit of the property to plot:')
prop_ul=input('Enter prop_ul limit of the property to plot:')
#prop= prop+4
title= ['Vp(km/s)','Vs(km/s)','Vp/Vs','Density(kg/m3)','wt%','vol%','mol%']

title_save= ['Vp','Vs','Vp_Vs','Density','Weight (%)','Volume (%)','mol%']
#mat_list =[90,97,99]
mat_list=input('Enter Composition files code (e.g 90,99 etc) in sequence separated by comma:' )
print mat_list
dist=input('Enter distance (in km and should be multiple of the resolution along the profile) along which to plot property : ')
post_data=np.loadtxt('post_processing_output.dat',usecols=(0,1,2,3,4,5,6))


### Olivine
d_Ol= find_mineral('Ol',mat_list,prop)
d=np.where(d_Ol[:,0]==dist)
y=d_Ol[d[:],1]
x=d_Ol[d[:],2]
t=d_Ol[d[:],3]
fig  = plt.figure(figsize=(5,6))
ax1 = fig.add_subplot(111)
ax1.plot(x[0,:],y[0,:],label='Olivine',color='limegreen')

d_Gt= find_mineral('Gt',mat_list,prop)
d=np.where(d_Gt[:,0]==dist)
y=d_Gt[d[:],1]
x=d_Gt[d[:],2]
t=d_Gt[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Garnet',color='darkorange')

d_Opx= find_mineral('Opx',mat_list,prop)
d=np.where(d_Opx[:,0]==dist)
y=d_Opx[d[:],1]
x=d_Opx[d[:],2]
t=d_Opx[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Orthopyroxene',color='lightgreen')

d_Cpx= find_mineral('Cpx',mat_list,prop)
d=np.where(d_Cpx[:,0]==dist)
y=d_Cpx[d[:],1]
x=d_Cpx[d[:],2]
t=d_Cpx[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Clinopyroxene',color='darkgreen')
'''
d_Wad= find_mineral('Wad',mat_list,prop)
d=np.where(d_Wad[:,0]==dist)
y=d_Wad[d[:],1]
x=d_Wad[d[:],2]
t=d_Wad[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Wadsleyite')

d_Ph= find_mineral('Ph',mat_list,prop)
d=np.where(d_Ph[:,0]==dist)
y=d_Ph[d[:],1]
x=d_Ph[d[:],2]
t=d_Ph[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Plagioclase')

d_Sp= find_mineral('Sp',mat_list,prop)
d=np.where(d_Sp[:,0]==dist)
y=d_Sp[d[:],1]
x=d_Sp[d[:],2]
t=d_Sp[d[:],3]
ax1.plot(x[0,:],y[0,:],label='Spinel')
'''
#ax1.axhline(y=100,linestyle='dotted',color='black',linewidth=1)
ax1.grid(True,linestyle='dotted',color='black',linewidth=0.1)
plt.legend( fancybox=True,shadow=False, framealpha=1,loc='best',fontsize=legend_font_size)

#ax1.invert_xaxis()
ax1.set_title('HP_AZ Distance =' + str(dist) + ' km',fontsize=axes_label_font_size, fontweight='bold')
ax1.set_ylabel('Depth (km)',fontsize=axes_label_font_size, fontweight='bold')
ax1.set_xlabel(str(title_save[prop-1]),fontsize=axes_label_font_size, fontweight='bold')
plt.ylim([120,400])
#plt.xlim([3.9,5.5])
plt.xlim([prop_lw,prop_ul])

plt.xticks(fontsize=axes_label_font_size)
plt.yticks(fontsize=axes_label_font_size)
#plt.ylim([0,100])
ax1.invert_yaxis()
#plt.tight_layout(pad=1.2, w_pad=1.2, h_pad=2.0)
s=(str(cwd+ '/' + str(title_save[prop-1]) + str(dist) + 'km_plot.png'))
savefig(s, dpi=400) 
#plt.tight_layout()

plt.show()
