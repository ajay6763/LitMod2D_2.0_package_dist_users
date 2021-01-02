#!/bin/sh
#####################################################################
#################	September 2019	@ajay	#####################
#####################################################################
##### This script calculated P-wave reciever functions          #####
##### using  CPS  (Herman 2013). CPS needs to installed and in 
####   path for     #####
##### this script to work.                                      #####
##### It produces RF at each distance using velocity CPS input  #####
##### file							#####
##### Inputs: 1. Gaussian width, 2. rap parameter              	##### 
##### provided.							#####
##### Minimum and Maximum period is 25 s  to 250 s              ##### 
#####################################################################
#################	Tasks	    @ajay	#####################
#####################################################################
##### Plotting part is independent				#####
#####  automate it -> @ajay6763					#####
#####################################################################
#####################################################################
echo 
echo "#####################################################################"
echo "#####################################################################"
echo " Now we will calculate synthetic P-wave reciver functions            "
echo
echo " Enter the ray parameter (usual value is 0.06) :" 
read P
#P=0.06
echo
echo " Enter the Gaussian width  (e.g. 2.5, 2.0) :" 
read ALP
#ALP=2.5
for file in *_vel.dat
do
hrftn96 -P -ALP ${ALP} -DT 0.1 -D 5. -RAYP ${P} -M $file  -2 hr	-NSAMP 1500	
mv hrftn96.sac ${file}_${ALP}.eqr
done



