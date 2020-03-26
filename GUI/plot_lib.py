#### ploting liberary
import matplotlib as plt
from matplotlib.patches import Polygon
import matplotlib as mpl
import numpy as np
import os
from matplotlib.widgets import MultiCursor
from matplotlib import style
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import *
import matplotlib
matplotlib.rcParams.update({'errorbar.capsize': 2})
origin = 'lower'
label_size = 4
marker_size =2
mpl.rcParams['xtick.labelsize'] = label_size 
#import smoothListGaussian
#style.use('classic') #### style o
def smoothListGaussian(list,degree=5):  

     window=degree*2-1  

     weight=np.array([1.0]*window)  

     weightGauss=[]  

     for i in range(window):  

         i=i-degree+1  

         frac=i/float(window)  

         gauss=1/(np.exp((4*(frac))**2))  

         weightGauss.append(gauss)  

     weight=np.array(weightGauss)*weight  

     smoothed=[0.0]*(len(list)-window)  

     for i in range(len(smoothed)):  

         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  

     return smoothed  

def plot_lib(x_nodes,l,b,topo_in_file,bouger_in_file,geoid_in_file,FA_in_file,heatflux_in_file,anomaly,c,l_nodes,b_nodes,mat,max_depth,
  label_size,marker_size,anomaly_x,anomaly_y,anomaly_amount,anomaly_type,anomaly_compo,anomaly_color):
  #
  # Making string list of the observables file input
  #
  #data_list=[];data_list.extend([topo_in_file,bouger_in_file,FA_in_file,geoid_in_file,heatflux_in_file])
  #
  ## Figuring out what input files are given
  #name_len=[ len(data_list[k]) for k in range(len(data_list)) ];

  #for i in range(len(name_len)):
  #  if name_len[i] != 0 :
      # plot it

  #  else
  #    pass
  
  try:
    topo=np.loadtxt(topo_in_file)#f=np.loadtxt(self.digitized,comments=">")
    max_depth = -max_depth
    min_depth = max(topo[:,1]/1000)
  except Exception:
    max_depth = -400
    min_depth = 3
  # Initializing the figure
  #
  fig  = plt.figure()
  #gs = gridspec.GridSpec(6, 1, width_ratios=[3, 1]) 
  gs = gridspec.GridSpec(6, 1,height_ratios=[1,1,1,1,1,3])
  #
  # ploting the bodies geometry
  #
  ax1 = plt.subplot(gs[5]) #fig.add_subplot(616)
  for i in range(len(l_nodes)):
                ax1.plot(l_nodes[i],b_nodes[i],color='black',lw=1,marker=None)
                poly = Polygon(list(zip(l_nodes[i],b_nodes[i])),facecolor=c[i],label=mat[i], lw=1)
                ax1.add_patch(poly)
  for i in range(len(anomaly_x)):
    ax1.plot(anomaly_x[i],anomaly_y[i],color="green",lw=1,marker='o',markersize=4)
    poly = Polygon(list(zip(anomaly_x[i],anomaly_y[i])),facecolor=anomaly_color[i],label="anomaly", lw=1)
    l_text=min(anomaly_x[i])+ (max(anomaly_x[i])-min(anomaly_x[i]))/2
    b_text=min(anomaly_y[i])+ (max(anomaly_y[i])-min(anomaly_y[i]))/2
    ax1.add_patch(poly)
    ax1.text(l_text, b_text, str(int(i+1))+" "+str(anomaly_type[i])+" Material"
      +str(anomaly_compo[i])+""+str(anomaly_amount[i]), fontsize=8,fontstyle='italic',color='black',bbox=dict(facecolor='white', alpha=0.9))
  #ax1.plot(l[i][:],b[i][:],color='grey',linewidth=0.5)
  ax1.grid(True)
  ax1.tick_params(labelsize=10)
  plt.ylim((max_depth,min_depth))
  plt.title('Profile Geometry')
  plt.xlabel('Distance ($km$)')
  plt.ylabel('Depth ($km$)')
  #
  # plotting the topography
  # 
  ax2 = plt.subplot(gs[4])#fig.add_subplot(611)    
  # topography data
  f_out=open("topo_out.dat","r")
  data_out=f_out.readlines()
  try:
    f_in=open(topo_in_file,"r") #f_in=open("test_top.dat","r")#f_in=open("test_top.dat","r")
    data_in=f_in.readlines()
    f_in.close()
  except Exception:
    data_in=data_out
  profile_out=[]
  profile_in=[]
  if anomaly==1:
    topo_in=[]
    topo_error=[]
    topo_out_coupled=[]
    topo_out_uncoupled=[]
    f_out.close()
    ## storing data into local variables
    for i  in range(len(data_in)):
      try:
        data_temp=data_out[i].split()
      except IndexError:
        data_temp=[]
        pass
      try:
        profile_out.append(float(str(data_temp[0])))
        topo_out_coupled.append(float(str(data_temp[1])))
        topo_out_uncoupled.append(float(str(data_temp[2])))
      except ValueError:
        data_temp=[]
        pass
      try:
        data_temp=data_in[i].split()
      except IndexError:
        data_temp=[]
        pass
      try:
        profile_in.append(float(str(data_temp[0])))
        topo_in.append(float(str(data_temp[1])))
        try:
          topo_error.append(float(str(data_temp[2])))
        except:
          topo_error.append(0)
      except ValueError:
        data_temp=[]
        pass    
    plt.hold(True)
    ax2.errorbar(profile_in,topo_in,yerr=topo_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')
    ax2.plot(profile_out,topo_out_coupled,'k-',label='Coupled' )
    ax2.plot(profile_out,topo_out_uncoupled, 'r-',label='Uncoupled')
    '''
    ax2_=ax2.twinx()
    diff=[]
    diff=np.subtract(topo_out_uncoupled,topo_in)
    ax2_.plot(profile_out,diff,color='g',label='Diff')
    ax2_.set_ylabel('Diff', color='g')
    ax2_.tick_params('y', colors='g')
    '''
    #ax2.plot(profile_in,topo_in, 'b-',label='Observed')
    ax2.grid(True)
    ax2.set_title('Elevation ')
    ax2.xaxis.set_ticklabels([])
    ax2.tick_params(labelsize=10)
    #plt.xlabel('Distance (Km)')
    plt.ylabel('Elevation ($m$)')
    plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
    #ax2.legend(handles, ['input','Computed-coupled','Computed-uncoupled'])
  else:
    topo_in=[]
    topo_error=[]
    topo_out=[]
    f_out.close()
    ## storing data into local variables
    for i  in range(len(data_in)):
      try:
        data_temp=data_out[i].split()
      except IndexError:
        pass
      try:
        profile_out.append(float(str(data_temp[0])))
        topo_out.append(float(str(data_temp[1])))
      except ValueError:
        passa
      try:
        data_temp=data_in[i].split()
      except IndexError:
        pass
      try:
        profile_in.append(float(str(data_temp[0])))
        topo_in.append(float(str(data_temp[1])))
        try:
          topo_error.append(float(str(data_temp[2])))
        except:
          topo_error.append(0)
      except ValueError:
        pass   
    #topo_misfit= (topo_in - topo_out)/len(topo_in)
    #print topo_misfit 
    plt.hold(True)
    #ax2.plot(profile_in,topo_in, 'b-',label='Observed')
    ax2.errorbar(profile_in,topo_in,yerr=topo_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')    
    ax2.plot(profile_out,topo_out,'r-',label='Calculated')
    '''
    ax2_=ax2.twinx()
    diff=[]
    diff=np.subtract(topo_out,topo_in)
    ax2_.plot(profile_out,diff,color='g',label='Diff')
    ax2_.set_ylabel('Diff', color='g')
    ax2_.tick_params('y', colors='g')
    '''
    #ax2.fill(profile_in,topo_in+topo_error,zorder=10)
    ax2.grid(True)
    ax2.set_title('Elevation')
    ax2.xaxis.set_ticklabels([])
    ax2.tick_params(labelsize=10)
    #plt.xlabel('Distance (Km)')
    plt.ylabel('Elevation ($m$)')
    plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
  #  
  # bouguer data
  #
  f_out=open("bouguer_out.dat","r")
  data_out=f_out.readlines()
  try:
    f_in=open(bouger_in_file,"r") #f_in=open("test_bou.dat","r")
    data_in=f_in.readlines()
    f_in.close()
  except Exception:
    data_in=data_out
  profile_out=[]
  profile_in=[]
  boug_in=[]
  boug_out=[]
  boug_error=[]
  for i  in range(len(data_in)):
    try:
      data_temp=data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      boug_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      boug_in.append(float(str(data_temp[1])))
      try:
          boug_error.append(float(str(data_temp[2])))
      except:
          boug_error.append(0.0)
    except ValueError:
      pass  

  ax3 = plt.subplot(gs[3])#fig.add_subplot(612)
  plt.hold(True)
  #ax3.plot(profile_in,boug_in, 'b-')
  ax3.errorbar(profile_in,boug_in,yerr=boug_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')
  ax3.plot(profile_out,boug_out,'r-',label='Calculated')
  '''
  ax3_=ax3.twinx()
  diff=[]
  diff=np.subtract(boug_out,boug_in[range(len(boug_out))])
  ax3_.plot(profile_out,diff,color='g',label='Diff')
  ax3_.set_ylabel('Diff', color='g')
  ax3_.tick_params('y', colors='g')
  '''
  ax3.grid(True)
  ax3.set_title('Bouguer')
  ax3.xaxis.set_ticklabels([])
  ax3.tick_params(labelsize=10)
  #plt.xlabel('Distance (Km)')
  plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
  plt.ylabel('Bouguer ($mGal$)')
  f_out.close()
  
  #
  # geoid data
  #
  f_out=open("geoid_out.dat","r")
  data_out=f_out.readlines()
  try:
    f_in=open(geoid_in_file,"r")#f_in=open("test_geo.dat","r")
    data_in=f_in.readlines()
    f_in.close()
  except Exception:
    data_in=data_out
  profile_out=[]
  profile_in=[]
  geoid_in=[]
  geoid_out=[]
  geoid_error=[]
  for i  in range(len(data_in)):
    try:
      data_temp=data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      geoid_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      geoid_in.append(float(str(data_temp[1])))
      try:
        geoid_error.append(float(str(data_temp[2])))
      except IndexError:
        geoid_error.append(0.0)        
    except ValueError:
      pass  
  #geoid_in=np.loadtxt(geoid_in_file)
  ax4 = plt.subplot(gs[2])#fig.add_subplot(614)
  plt.hold(True)
  ax4.plot(profile_in,geoid_in, 'b-')
  ax4.errorbar(profile_in,geoid_in,yerr=geoid_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')
  '''
  ax4_=ax4.twinx()
  diff=[]
  diff=np.subtract(geoid_out,geoid_in[range(len(geoid_out))])
  ax4_.plot(profile_out,diff,color='g',label='Diff')
  ax4_.set_ylabel('Diff', color='g')
  ax4_.tick_params('y', colors='g')
  '''
  #test=[]
  #test=smoothListGaussian(geoid_in,degree=1)
  #ax4.plot(profile_in[0:len(test)],test,fmt='o-',markersize=marker_size,color='black',label='Observed--9')
  '''
  try:
    ax4.errorbar(geoid_in[:,0],geoid_in[:,1],yerr=geoid_in[:,2],fmt='o-',markersize=marker_size,ecolor='b',label='8')
    plt.hold(True)
  except Exception:
    pass
  try:
    ax4.errorbar(geoid_in[:,0],geoid_in[:,3],yerr=geoid_in[:,4],fmt='o-',markersize=marker_size,ecolor='k',label='9')
    plt.hold(True)
  except Exception:
    pass
  try:
    ax4.errorbar(geoid_in[:,0],geoid_in[:,5],yerr=geoid_in[:,6],fmt='o-',markersize=marker_size,ecolor='m',label='10')
    plt.hold(True)
  except Exception:
    pass
  '''
  plt.hold(True)
  ax4.plot(profile_out,geoid_out,'r-', label='Calculated')
  ax4.grid(True)
  ax4.set_title('Geoid')
  ax4.xaxis.set_ticklabels([])
  ax4.tick_params(labelsize=10)
  plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
  #plt.xlabel('Distance (Km)')
  plt.ylabel('Geoid ($m$)')
  f_out.close()
  #
  # heat flux data
  #
  f_out=open("SHF_out.dat","r")
  data_out=f_out.readlines()
  f_out.close()
  try:
    f_in=open(heatflux_in_file,"r")#f_in=open("fluxout.dat","r")
    data_in=f_in.readlines()
    f_in.close()
  except Exception:
    data_in=data_out
  profile_out=[]
  profile_in=[]
  flux_in=[]
  flux_out=[]
  flux_error=[]
  for i  in range(len(data_out)):
    try:
      data_temp=data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      flux_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass
  for i  in range(len(data_in)):
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      flux_in.append(float(str(data_temp[1])))
      try:
        flux_error.append(float(str(data_temp[2])))
      except IndexError:
        flux_error.append(0.0)
    except ValueError:
      pass
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass  
  

  ax5 = plt.subplot(gs[0])#fig.add_subplot(615)
  plt.hold(True)
  ax5.errorbar(profile_in,flux_in,yerr=flux_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')
  ax5.plot(profile_out,flux_out,'r-',label='Calculated')
  ax5.grid(True)
  ax5.set_title('Heat Flux')
  ax5.xaxis.set_ticklabels([])
  ax5.tick_params(labelsize=10)
  #plt.xlabel('Distance (Km)')
  plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
  plt.ylabel('Heat-flux ($mW/m^2$)')

  #plt.tight_layout()

  #plt.tight_layout()  
  #
  # FA data
  #
  f_out=open("FA_out.dat","r")
  data_out=f_out.readlines()
  f_out.close()
  try:
    f_in=open(FA_in_file,"r")#f_in=open("fluxout.dat","r")
    data_in=f_in.readlines()
    f_in.close()
  except Exception:
    data_in=data_out
  profile_out=[]
  profile_in=[]
  FA_in=[]
  FA_out=[]
  FA_error=[]
  for i  in range(len(data_in)):
    try:
      data_temp=data_out[i].split()
    except IndexError:
      pass
    try:
      profile_out.append(float(str(data_temp[0])))
      FA_out.append(float(str(data_temp[1])))
    except ValueError:
      pass
    try:
      data_temp=data_in[i].split()
    except IndexError:
      pass
    try:
      profile_in.append(float(str(data_temp[0])))
      FA_in.append(float(str(data_temp[1])))
      try:
        FA_error.append(float(str(data_temp[2])))
      except IndexError:
        FA_error.append(0.0)
    except ValueError:
      pass  
  ax6 = plt.subplot(gs[1])#fig.add_subplot(613)
  plt.hold(True)
  #ax6.plot(profile_in,FA_in, 'b-')
  ax6.errorbar(profile_in,FA_in,yerr=FA_error,fmt='o-',markersize=marker_size,ecolor='b',label='Observed')
  ax6.plot(profile_out,FA_out,'r-',label='Calculated')
  '''
  ax6_=ax6.twinx()
  diff=[]
  diff=np.subtract(FA_out,FA_in[range(len(FA_out))])
  ax6_.plot(profile_out,diff,color='g',label='Diff')
  ax6_.set_ylabel('Diff', color='g')
  ax6_.tick_params('y', colors='g')
  '''
  ax6.grid(True)
  ax6.set_title('Free-Air')
  ax6.xaxis.set_ticklabels([])
  ax6.tick_params(labelsize=10)
  plt.legend( fancybox=True, framealpha=0.5,fontsize=8,loc='best')
  #plt.xlabel('Distance (Km)')
  plt.ylabel('FA (mGal)')
  #plt.tight_layout()

  multi = MultiCursor(fig.canvas,(ax1,ax2,ax3,ax4,ax5,ax6), color='r', lw=1)


  ##############33 temperature profile plotting
  gs_ = gridspec.GridSpec(4, 1)
  #fig_out = plt.figure()
  fig_out  = plt.figure()
  #ax_temp = fig_out.add_subplot(411)
  ax_temp =plt.subplot(gs_[0]) #fig_out.add_subplot(411)
  f=open("tempout.dat","r")
  data_out=f.readlines()
  temp_x=[]
  temp_y=[]
  temp_z=[]
  for i  in range(2,len(data_out)):
    data_temp=data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z.append(float(str(data_temp[2])))
  xlist=[]
  ylist = np.array(temp_y[0:96])
  z =temp_z #np.array(temp_z)
  Z=[]
  temp=0
  for i in range(0,x_nodes):
    tmp = []
    xlist.append(temp_x[temp])
    for j in range(0,96):
      tmp.append(z[temp])
      temp=temp+1
    Z.append(tmp)
  X,Y = np.meshgrid(ylist, xlist)   
  levels = np.linspace(0, 50, 1600)
  
  
  for i in range(len(l)):
    ax_temp.plot(l[i][:],b[i][:],color='grey',linewidth=2)
  levels =[300,500,600,1000,1200,1300,1450,1550,1650]
  CS3 = plt.contourf(Y, X, Z, levels,extend='both',cmap=plt.cm.get_cmap('coolwarm',20))
  CS3.cmap.set_under('white')
  CS3.cmap.set_over('black')
  #levels = np.arange(min(temp_z), max(temp_z)+100, 100)
  CS4 = plt.contour(Y, X, Z, levels,linewidths=0.2,linestyles='dashed',colors='black')
  plt.clabel(CS4,inline=5,inline_spacing=10,fontsize=10,fontstyles='arial',fontweight='bold',fmt='%1.0f',colors='black')
  plt.colorbar(CS3,label='$^oC$')
  plt.title('Temperature',fontsize=10, fontweight='bold')
  """
  cp = plt.contourf(Y,X,Z)
  plt.colorbar(cp)
  #cpc=plt.contour(Y,X,Z,levels) # predefinrd number of contour
  CS = plt.contour(Y, X, Z, 30 ,linewidths=0.5,colors='grey') # automatice number of contours
  plt.clabel(CS,inline=1,inline_spacing=0,fontsize=12,fmt='%1.0f',colors='black')
  #plt.clabel(CS, inline=1, fontsize=10,color="black") # contour line labels
  """
  #plt.title('Temperature profile')
  #plt.xlabel('Distance ($km$)')
  plt.tick_params(labelsize=10)
  plt.ylabel('Depth ($km$)')

  ##############33 density profile plotting
  f=open("dens_node2.dat","r")
  data_out=f.readlines()
  temp_x=[]
  temp_y=[]
  temp_z=[]
  for i  in range(len(data_out)):
    data_temp=data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z.append(float(str(data_temp[2])))
  xlist=[]
  ylist = np.array(temp_y[0:95])
  z =temp_z #np.array(temp_z)
  Z=[]
  temp=0
  for i in range(0,x_nodes):
    tmp = []
    xlist.append(temp_x[temp])
    for j in range(0,95):
      tmp.append(z[temp])
      temp=temp+1
    Z.append(tmp)
  X,Y = np.meshgrid(ylist, xlist)   
  levels = np.linspace(0, 50, 1600)
  
  #ax_dens = fig_out.add_subplot(412)
  ax_dens =plt.subplot(gs_[1])
  for i in range(len(l)):
    ax_dens.plot(l[i][:],b[i][:],color='grey',linewidth=2)
  levels = np.arange(2800, 3800, 50)
  #levels = np.arange(3000, 3800, 20)
  CS3 = plt.contourf(Y, X, Z, levels,extend='both',cmap=plt.cm.get_cmap('rainbow',20),origin=origin)
  CS3.cmap.set_under('white')
  CS3.cmap.set_over('black')
  levels = np.arange(3000, 3800, 50)
  CS4 = plt.contour(Y, X, Z,levels,linewidths=0.2,linestyles='dashed',colors='grey')
  plt.clabel(CS4,inline=True,inline_spacing=2,fontsize=10,fontstyles='arial',fontweight='normal',fmt='%1.0f',colors='black')
  plt.colorbar(CS3,label='$kg/m^3$')
  plt.title('Density',fontsize=10, fontweight='bold')
  plt.tick_params(labelsize=10)
  """  
  cp = plt.contourf(Y,X,Z)
  plt.colorbar(cp)3
  #cpc=plt.contour(Y,X,Z,levels) # predefinrd number of contour
  CS = plt.contour(Y, X, Z, 30 ,linewidths=0.5,colors='grey') # automatice number of contours
  plt.clabel(CS,inline=1,inline_spacing=0,fontsize=12,fmt='%1.0f',colors='black')
  #plt.clabel(CS, inline=1, fontsize=10,color="black") # contour line labels
  """
  #plt.title('Density profile')
  #plt.xlabel('Distance ($km$)')
  plt.ylabel('Depth ($km$)')
 

  ##############33 Vp profile plotting
  #f=open("velatten.dat","r")
  f=open("velocities.dat","r")
  data_out=f.readlines()
  temp_x=[]
  temp_y=[]
  temp_z_vp=[]
  temp_z_vs=[]
  for i  in range(len(data_out)):
    data_temp=data_out[i].split()
    temp_x.append(float(str(data_temp[0])))
    temp_y.append(float(str(data_temp[1])))
    temp_z_vp.append(float(str(data_temp[2])))
    temp_z_vs.append(float(str(data_temp[3])))
  xlist=[]
  ylist = np.array(temp_y[0:95])
  z_vp =temp_z_vp #np.array(temp_z)
  z_vs =temp_z_vs
  Z_vp=[]
  Z_vs=[]
  temp=0
  for i in range(0,x_nodes):
    tmp1 = []
    tmp2 = []
    xlist.append(temp_x[temp+1])
    for j in range(0,95):
      tmp1.append(z_vp[temp])
      tmp2.append(z_vs[temp])
      temp=temp+1
    Z_vp.append(tmp1)
    Z_vs.append(tmp2)
  X,Y = np.meshgrid(ylist, xlist)   
  levels = np.linspace(0, 50, 1600)
  
  #ax_vp = fig_out.add_subplot(413)
  ax_vp =plt.subplot(gs_[2])
  for i in range(len(l)):
    ax_vp.plot(l[i][:],b[i][:],color='grey',linewidth=2)
  levels = np.arange(7.2, 8.8, 0.01)
  CS3 = plt.contourf(Y, X, Z_vp, levels,extend='both',cmap=plt.cm.get_cmap('RdYlBu',20),origin=origin)
  #CS3 = plt.contourf(Y, X, Z_vp, levels,extend='both',cmap=cmap,origin=origin)
  CS3.cmap.set_under('white')
  CS3.cmap.set_over('black')
  levels = np.arange(7.7, 8.8, 0.05)
  CS4 = plt.contour(Y, X, Z_vp, levels,linewidths=0.2,linestyles='dashed',colors='black')
  plt.clabel(CS4,inline=5,inline_spacing=10,fontsize=10,fmt='%1.1f',colors='black')
  plt.colorbar(CS3,label='$km/s$')
  plt.tick_params(labelsize=10)
  plt.title('P wave velocity',fontsize=10, fontweight='bold')
  #plt.xlabel('Distance ($km$)')
  plt.ylabel('Depth ($km$)') 


  #ax_vs = fig_out.add_subplot(414)
  ax_vs =plt.subplot(gs_[3])


  for i in range(len(l)):
    ax_vs.plot(l[i][:],b[i][:],color='grey',linewidth=2)
  levels = np.arange(4.2, 4.8, 0.01)
  #CS3 = plt.contourf(Y, X, Z_vs, levels,extend='both')
  CS3 = plt.contourf(Y, X, Z_vs, levels,extend='both',cmap=plt.cm.get_cmap('RdYlBu',20),origin=origin)
  CS3.cmap.set_under('white')
  CS3.cmap.set_over('black')
  #levels = np.arange(4.1, 4.8, 0.05)
  levels = np.arange(4.3, 5.1, 0.05)
  CS4 = plt.contour(Y, X, Z_vs, levels,linewidths=0.2,linestyles='dashed',colors='black')
  plt.clabel(CS4,inline=5,inline_spacing=10,fontsize=10,fmt='%1.1f',colors='black')
  plt.colorbar(CS3,label='$km/s$')
  plt.title('S wave velocity',fontsize=10, fontweight='bold')
  plt.xlabel('Distance ($km$)')
  plt.ylabel('Depth ($km$)')
  plt.tick_params(labelsize=10)
  #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
  #plt.tight_layout()
  savefig('Temp_Dens_velocities.png', dpi=300) 
  

  ###########################3
  ### calculating the misfit
  
  
  #boug_misfit= (bouger_in - bouguer_out)/len(bouger_in)
  #print boug_misfit
  #geoid_misfit=(geoid_in - geoid_out )/len(geoid_in)
  #print geoid_misfit

  ##################
  plt.show()