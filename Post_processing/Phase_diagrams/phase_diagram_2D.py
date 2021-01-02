#### ploting liberary
import matplotlib as plt
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches
import numpy as np
import os
from matplotlib.widgets import MultiCursor
from matplotlib import style
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import *
import math
cwd = os.getcwd()

import matplotlib

matplotlib.rcParams.update({'font.size': 16})

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
#######################
## functions to find the index
'''
def lookup(T_LitMod,P_LitMod,T,P):
	index=[]
#	dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
	dist=np.array((T-T_LitMod)**2-(P-P_LitMod)**2);
	index=dist.argmin(); 
	#print index
	return index
'''
def lookup(T_LitMod,P_LitMod,T,P):
	index=[]
#	dist=np.array((T[:]-T_LitMod)**2-( P[:]-P_LitMod)**2)
	dist=np.array(sqrt((T-T_LitMod)**2+(P-P_LitMod)**2));
	index=dist.argmin(); 
	#print index, T_LitMod,P_LitMod
	return index	
def find_mineral(mineral,mat_list,d):
	Ol=[]
	X=[]
	Z=[]
	T_out=[]
	d=d+1
	for i in range(len(mat_list)):
		s=str(str(mat_list[i])+"_profile")
		data=np.loadtxt(s)
		#T P Vp Vs Vp/Vs density wt% vol%
		s=str(str( mineral +"_"+ str(mat_list[i]))+"_FULL.txt")
		Ol_data=np.loadtxt(s)
		if len(Ol_data) ==0:
			pass
		else:
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
prop=input('Enter what property you want to plot: \n 1- P-wave velocities \n 2- S-wave velocities \n 3- Vp/Vs \n 4- Density \n 5- Weight % \n 6- Volume % \n 7- mole % \n  ')

title= ['Vp(km/s)','Vs(km/s)','Vp/Vs','Density(kg/m3)','wt%','vol%','mole%']

title_save= ['Vp','Vs','Vp_Vs','Density','wt%','vol%','mole%']
mat_list =input('Enter Composition files code (e.g 90,99 etc) in sequence separated by comma:' )

d=[]
d= find_mineral('Ol',mat_list,prop)
#d.sort(axis=-2)
#d[np.lexsort(np.fliplr(d).T)]
f=open("Olivine.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()
fig  = plt.figure(figsize=(25,14))
#fig2  = plt.figure(figsize=(25,14))

ax1 = fig.add_subplot(421)
#ax1_ = fig2.add_subplot(321)
#ax1 = plt.subplot(321)
size=np.divide(d[:,2],max(d[:,2]))
size=80
symbol='s'
symbol_edge_color='None'
symbol_edge_line=0.1
cm = plt.cm.get_cmap('Purples')
sc = ax1.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax1.set_title('Olivine')
ax1.set_xlabel('Distance (km)')
ax1.set_ylabel('Depth (km)')
ax1.invert_yaxis()
'''
cm_ = plt.cm.get_cmap('Purples')
sc_ = ax1_.scatter(d[:,0], d[:,1], c=d[:,3], vmin=min(d[:,3]), vmax=max(d[:,3]), s=size, cmap='Spectral'_,marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar_=plt.colorbar(sc_)
cbar_.ax.set_title(title[prop-1])
ax1_.set_title('Olivine')
ax1_.set_xlabel('Distance (km)')
ax1_.set_ylabel('Depth (km)')
ax1_.invert_yaxis()
'''



############################################
############ Gt
d=[]
d= find_mineral('Gt',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("Garnet.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()

#ax2 = plt.subplot(322)
ax2 = fig.add_subplot(422)
#size=np.divide(Ol*10,max(Ol))
sc = ax2.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax2.set_title('Garnet')
ax2.set_xlabel('Distance (km)')
ax2.set_ylabel('Depth (km)')
ax2.invert_yaxis()
#################################################
############ Opx

d=[]
d= find_mineral('Opx',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("Orthopyroxene.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()

#ax3 = plt.subplot(323)
ax3 = fig.add_subplot(423)
#size=np.divide(Ol*10,max(Ol))
sc = ax3.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax3.set_title('Orthopyroxene')
ax3.set_xlabel('Distance (km)')
ax3.set_ylabel('Depth (km)')
ax3.invert_yaxis()

################################################3
###############3 CpX

d=[]
d= find_mineral('Cpx',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("Clinopyroxene.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()

#size=np.divide(Ol*10,max(Ol))
#ax4 = plt.subplot(324)
ax4 = fig.add_subplot(424)
sc = ax4.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax4.set_title('Clinopyroxene')
ax4.set_xlabel('Distance (km)')
ax4.set_ylabel('Depth (km)')
ax4.invert_yaxis()

##################3
### pyroxene
d=[]
d= find_mineral('C2c',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("pyroxene.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()
#size=np.divide(Ol*10,max(Ol))
#ax5 = plt.subplot(325)
ax5 = fig.add_subplot(425)
sc = ax5.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax5.set_title('Pyroxene')
ax5.set_xlabel('Distance (km)')
ax5.set_ylabel('Depth (km)')
ax5.invert_yaxis()



d=[]
d= find_mineral('Ph',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("Plagioclase.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()
#size=np.divide(Ol*10,max(Ol))
#ax5 = plt.subplot(325)
ax6 = fig.add_subplot(426)
sc = ax6.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax6.set_title('Plagioclase')
ax6.set_xlabel('Distance (km)')
ax6.set_ylabel('Depth (km)')
ax6.invert_yaxis()


d=[]
d= find_mineral('Sp',mat_list,prop)
#d.sort(axis=-2)
d[np.lexsort(np.fliplr(d).T)]
f=open("Spinel.dat","w")
for i in range(len(d)):
	f.writelines(" %f  %f  %f  \n " % (d[i,0],d[i,1],d[i,2]))
f.close()

#size=np.divide(Ol*10,max(Ol))
#ax6 = plt.subplot(326)
ax7 = fig.add_subplot(427)
sc = ax7.scatter(d[:,0], d[:,1], c=d[:,2],vmin=min(d[:,2]), vmax=max(d[:,2]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title(title[prop-1])
ax7.set_title('Spinel')
ax7.set_xlabel('Distance (km)')
ax7.set_ylabel('Depth (km)')
ax7.invert_yaxis()

'''
post_data=np.loadtxt('post_processing_output_crust.dat',usecols=(0,1,4,5))
vp_vs=np.divide(post_data[:,2],post_data[:,3],out=np.zeros_like(post_data[:,2]), where=post_data[:,3]!=0)
ax8 = fig.add_subplot(428)
sc = ax8.scatter(post_data[:,0], -post_data[:,1], c=vp_vs[:], vmin=min(vp_vs[:]), vmax=max(vp_vs[:]), s=size, cmap='Spectral',marker=symbol,edgecolors=symbol_edge_color,linewidth=symbol_edge_line)
cbar=plt.colorbar(sc)
cbar.ax.set_title('Vp/Vs')
ax8.set_title('Vp/Vs')
ax8.set_xlabel('Distance (km)')
ax8.set_ylabel('Depth (km)')
ax8.invert_yaxis()
'''


#plt.tight_layout(pad=0.2, w_pad=0.2)
plt.tight_layout(pad=1.2, w_pad=1.2, h_pad=2.0)
s=(str(cwd+ '/' + str(title_save[prop-1]) + '_plot.jpg'))
savefig(s, dpi=400) 
#plt.tight_layout()
plt.show()

