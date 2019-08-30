import numpy as np
length=input('Enter Profile length:  ')
reso=input('Enter resolution along length: ')
data=np.loadtxt("data_phase_vel",)
x_nodes=length/reso
prog='SURF96'
control=1
control = input("What waves you want Love ( 0 ) or Rayleigh ( 1 - default value )  : ")
if control==1:
	surf_wave='R'
elif control==0:
	surf_wave='L'
else:
	surf_wave='R'




control=1
control=input("Do you want phase ( 0 ) or group ( 1 - default ) velocities : ")
if control==1:
	vel_type='U'
elif control==0:
	vel_type='C'
else:
	vel_type='U'



mod=0
xx='T'
error= 0.0
time_period=[4 ,6, 8, 10, 12, 14 ,16 ,18 ,20 ,25 ,30, 35, 40, 45, 50, 59, 67, 77, 87, 100, 111, 125, 143,167]

#try:
for i in range(x_nodes):
	if data[i,0] == 'NaN':
		print 'yes'
		pass 
	elif data[i,0] > 0:
		#s=str(str(i)+"_"+str(surf_wave)+"_"+str(vel_type)+"_SURF96.inp")
		s=str(str(i)+"_SURF96.inp")
		f=open(s,"w")
		for j in range(len(time_period)):
			f.writelines("%s \t %s \t %s \t %s \t %d \t %8.4f \t %f \t %f \n" %(prog,surf_wave,vel_type,xx,mod,time_period[j], data[i,j],error))
		f.close()
	else:
		continue
#except Exception:
#	print("Something is wronge with data.\n Check if you really have dispersion data at each node in your profile. \n Input file for dispersion calculations are made upto the point where you have dispersion data.")
