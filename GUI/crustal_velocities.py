import numpy as np 
length=input('Enter Profile length:  ')
reso=input('Enter resolution along length:  ')
data=np.loadtxt('post_processing_output.dat',usecols=(0, 1,2,3,4,5,6))
vp_vs=1.73
### data format
#  X	Z	T 	P 	Vp 	Vs 	density 	Material mk??
## caluclating Vp in crust according to Brocher 2005
f=0
for i in range(len(data)):
	if data[i,5]==0:
		dens=data[i,6]/1000
		data[i,4]= 39.128*dens - 63.064*dens**2 + 37.083*dens**3 - 9.1819*dens**4 + 0.8215*dens**5
		#f=vp_vs -  data[i,1]/1000
		data[i,5]=data[i,4]/vp_vs
		### Accroding to equation 6 Brocher BSSA 2005
		#data[i,5]= 0.7858 -1.2344*data[i,4] + 0.7949*data[i,4]**2 - 0.1238*data[i,4]**3 + 0.0064*data[i,4]**4
		
		### Accroding to equation 7 Brocher BSSA 2005
		#data[i,5]= (data[i,4]-1.36)/1.16 # does not work

		### Accroding to equation 8 Brocher BSSA 2005
		#data[i,5]= 2.88 + 0.52*(data[i,4]-5.25)  # works but not perfect

	else:
		pass


f=open("velocities_crust_mantle.dat","w")
for i in range(len(data)):
	if i== len(data)-1:
		f.writelines(" %f  %f  %f  %f  %f  " % (data[i,0],-data[i,1],data[i,6]/1000, data[i,4],data[i,5],))
	else:
		f.writelines(" %f  %f  %f  %f  %f  \n " % (data[i,0],-data[i,1],data[i,6]/1000, data[i,4],data[i,5],)) 
f.close()
data=np.loadtxt("velocities_crust_mantle.dat",)
x_nodes=length/reso

for i in range(x_nodes):
	s=str(str(i)+"_vel.dat")
	f=open(s,"w")
	f.write("MODEL\n\
TEST MODEL \n\
ISOTROPIC \n\
KGS \n\
FLAT EARTH \n\
1-D \n\
CONSTANT VELOCITY \n\
LINE08 \n\
LINE09 \n\
LINE10 \n\
LINE11 \n \
HR         VP       VS         RHO    QP    QS    ETAP  ETAS FREFP FREFS \n")
	for j in range(i+ (94*i),i+ (94*i) + 94):
		#print  i,j				
			f.writelines("%f    \t  %f    \t  %f   \t  %f   \t   %s   \t  %s  \t %s \t  %s \t  %s  \t %s \n"
    	 	 %(data[j+1,1]-data[j,1] ,  data[j,3],   data[j,4],  data[j,2], "200.0",  "100", "0.0",  "0.0",  "1.0", "1.0"))
	#f.writelines("%f    \t  %f    \t  %f   \t  %f   \t   %s   \t  %s  \t %s \t  %s \t  %s  \t %s \n"
    # 	 %(data[-2,1]-data[-3,1] ,  data[-2,3],   data[-2,4],  data[-2,2], "200.0",  "100", "0.0",  "0.0",  "1.0", "1.0"))
	f.writelines("%f    \t  %f    \t  %f   \t  %f   \t   %s   \t  %s  \t %s \t  %s \t  %s  \t %s \n"
     	 %(data[-1,1]-data[-2,1] ,  data[-1,3],   data[-1,4],  data[-1,2], "200.0",  "100", "0.0",  "0.0",  "1.0", "1.0"))
	
	f.close()

	