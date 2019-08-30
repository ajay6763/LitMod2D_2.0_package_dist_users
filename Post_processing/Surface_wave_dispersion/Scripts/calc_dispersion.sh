#!/bin/bash
source=/home/akumar/owncloud/PHD/LITMOD_package_Linux_dist/Surface_wave_dispersion/Scripts
cp $source/modl.d ./
cp $source/sobs.d ./
cp $source/disp.d ./
for file  in *SURF96.inp
do
i=`echo $file| awk -F_ '{print $1}'`
cp $i"_vel.dat" modl.d
cp $i"_SURF96.inp" disp.d
surf96 39 
surf96 32 1 
surf96 36 1 
surf96 30 1 
surf96 17
surf96 1 
surf96 27 disp.out
mv disp.out $i"_SURF96.out"
done
