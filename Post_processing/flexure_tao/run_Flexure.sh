#!/bin/bash
echo "Enter length of the profile"
cp litmod.inp litmod_bu.inp
sed -e "24c   7    7    0    1" ./litmod_bu.inp > litmod.inp
/home/akumar/owncloud/PHD/LITMOD_package_Linux_dist/LITMOD_V4_LINUX
/home/akumar/owncloud/PHD/LITMOD_package_Linux_dist/flexure_tao/Pressure_Flexure.job
cp litmod_bu.inp litmod.inp
