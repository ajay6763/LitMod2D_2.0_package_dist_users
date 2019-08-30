import numpy as np 
### data format
#  X	Z	T 	P 	Vp 	Vs 	density 	Material mk??
## caluclating Vp in crust according to Brocher 2005
data=np.loadtxt("ak135f.txt",skiprows=0+1+2)

s=str("modl.d")
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


for j in range(len(data)-1):
	#print  i,j
	f.writelines("%s %f    \t  %f    \t  %f   \t  %f   \t   %f   \t  %f  \t %s \t  %s \t  %s  \t %s \n"
    	 %(str("f.writelines("),data[j+1,0]-data[j,0] ,  data[j,2],   data[j,3],  data[j,1], data[j,4],  data[j,5], "0.0",  "0.0",  "1.0", "1.0"))
f.close()

	