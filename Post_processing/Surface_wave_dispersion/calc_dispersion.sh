#!/bin/bash
#####################################################################
#################	September 2019	@ajay	#####################
#####################################################################
##### This script calculated dispersion curves along the profile
##### using  CPS  (Herman 2013). CPS needs to installed and in 
####   path for     #####
##### this script to work.                                      #####
##### It produces RF and dispersion curves at each distance     #####
##### point along the profile.					#####
##### Option for Love or Rayleigh Phase or Group velocities are ##### 
##### provided.							#####
##### Minimum and Maximum period is 25 s  to 250 s              ##### 
#####################################################################
#################	Tasks	    @ajay	#####################
#####################################################################
##### Plotting part is independent				#####
#####  automate it -> @ajay6763					#####
#####################################################################
#####################################################################
#source=$HOME/owncloud/PHD/LitMod_Methodology/LitMod2D_all_development/LitMod2D_2.0_package_dist_users/Post_processing
#$LitModHOME/Post_processing/Surface_wave_dispersion

####################################################################
###copy default disperions and model files
cp $LitModHOME/Post_processing/Surface_wave_dispersion/modl.d ./
cp $LitModHOME/Post_processing/Surface_wave_dispersion/sobs.d ./
cp $LitModHOME/Post_processing/Surface_wave_dispersion/disp.d ./


#####################################################################
### Asking for Love or Rayleigh Wave
ans=1
echo
echo 'What kind of sufrace wave you want (Rayleigh---R (default case) , Love--L):'
read surf_type
echo
echo 'Do you want Phase (C) or Group (U) velocities (default case is Group velocities):'
read disp_type
#####################################################################
### if Love is chosen than calculates Love wave dispersion both 
### Group (U in the name) and Phase (C in the name) 
#####################################################################

### Now for each distance node along the profile taking the velocity
### input file formatted for the CPS

for file  in *_vel.dat
do

#####################################################################
### getting the file identifier (distance along the profile)
### and copying the the model to modl.d (will be read by CPS)

i=`echo $file| awk -F_ '{print $1}'`
cp $i"_vel.dat" modl.d

#####################################################################
### Checking if Observed dispersion exists. If yes than using it 

if [  -f ${i}_${surf_type}_${disp_type}_SURF96.inp ] #|| [ -f  "${i}_L_C_SURF96.inp"  ]
then
    cp ${i}_${surf_type}_${disp_type}_SURF96.inp disp.d
else
#####################################################################
### If Observed dispersion doesn not exist than using dispersion curve
### for PREM model as observed (from GDM52 Ekstrome) 
    cp $LitModHOME/Post_processing/Surface_wave_dispersion/disp_${surf_type}_${disp_type}.d ./disp.d    
fi
#####################################################################
### Now calculating the dispersion curve
surf96 39 
surf96 32 1 
surf96 36 1 
surf96 30 1 
surf96 17
surf96 1 
surf96 27 disp.out
#####################################################################
### making the syntheic dispersion file
### format: {distance}km_{wave type}_{Phase/Group}
mv disp.out ${i}_${surf_type}_${disp_type}_SURF96.out 
done
