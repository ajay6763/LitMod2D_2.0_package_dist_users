#!/bin/bash
#echo "Enter length of the profile"
path=`echo $LitModHOME`
cp litmod.inp litmod_bu.inp
total_bodies=`awk 'FNR==24  {print $1}' litmod.inp`
sed -e "24c   ${total_bodies}    ${total_bodies}    0    1" ./litmod_bu.inp > litmod.inp
/home/temporary/owncloud/PHD/LitMod_Methodology/LitMod2D_all_development/LitMod2D_2.0_package_dist_users/GUI/LITMOD_V4_LINUX
/home/temporary/owncloud/PHD/LitMod_Methodology/LitMod2D_all_development/LitMod2D_2.0_package_dist_users/Post_processing/flexure_tao/Pressure_Flexure.job
cp litmod_bu.inp litmod.inp
#/home/temporary/owncloud/PHD/LitMod_Methodology/LitMod2D_all_development/LitMod2D_2.0_package_dist_users/GUI/LITMOD_V4_LINUX
