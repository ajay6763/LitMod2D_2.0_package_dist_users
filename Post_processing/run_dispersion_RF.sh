#!/bin/bash
#####################################################################
#################	September 2019	@ajay	#####################
#####################################################################
##### This script uses python script ot calculate velocities    #####
##### in the crust and produces velocity input files for CPS    #####
##### (Herman 2013). CPS needs to installed and in path for     #####
##### this script to work.                                      #####
##### It produces RF and dispersion curves at each distance     #####
##### point along the profile.					#####


#####################################################################
#################	Tasks	    @ajay	#####################
#####################################################################
##### Plotting part is independent				#####
#####  automate it -> @ajay6763					#####
#####################################################################
#####################################################################
python $LitModHOME/Post_processing/crustal_velocities_full.py

if [ -d "./Surface_Wave_Dispersion_Curves" ]; then
  ### Take action if $DIR exists ###
  echo "You have a Surface_Wave_Dispersion_Curves directory. It will be overwritten."
else
  ###  Control will jump here if $DIR does NOT exists ###
  mkdir ./Surface_Wave_Dispersion_Curves
fi


cp *_vel.dat ./Surface_Wave_Dispersion_Curves
cd ./Surface_Wave_Dispersion_Curves
${LitModHOME}/Post_processing/Surface_wave_dispersion/calc_dispersion.sh
cd ./../



if [ -d "./Receiver_Functions" ]; then
  ### Take action if $DIR exists ###
  echo "You have a Receiver_Functions directory. It will be overwritten."
else
  ###  Control will jump here if $DIR does NOT exists ###
  mkdir ./Receiver_Functions
fi

cp *_vel.dat ./Receiver_Functions
cd ./Receiver_Functions
${LitModHOME}/Post_processing/RF/Synthetic_RF/syn_P_RF.sh
for file in *.eqr
do
sac2xy $file $file.xyz
done
cd ./../

rm *_vel.dat

#/home/akumar/owncloud/PHD/LITMOD_package_Linux_dist/Surface_wave_dispersion/Scripts/disp_plot_along_profile.gmt


