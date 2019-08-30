#!/bin/bash
source=/home/akumar/owncloud/PHD/LITMOD_package_Linux_dist_users/Post_processing/Surface_wave_dispersion
cp $source/modl.d ./
cp $source/sobs.d ./
cp $source/disp.d ./
for file  in *_vel.dat
do
i=`echo $file| awk -F_ '{print $1}'`
cp $i"_vel.dat" modl.d
if [  -f "${i}_SURF96.inp" ]
then
    cp $i"_SURF96.inp" disp.d
else
    cp $source/disp.d ./    
fi
surf96 39 
surf96 32 1 
surf96 36 1 
surf96 30 1 
surf96 17
surf96 1 
surf96 27 disp.out
mv disp.out $i"_SURF96.out"
done
