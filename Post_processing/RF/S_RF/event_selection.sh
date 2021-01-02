#!/bin/bash
ub=50;
lb=30;
for dir in [12]*
do
cd $dir
dist_=`sachead c45_$dir'_G_HYB.BHZ' GCARC`;
dist=`echo $dist_|cut -f1 -d"."`;
if [ "$dist" -le "50" ];then
cp -r ./../$dir /DATABACKUP/STUDENT_DATA/ajay_data/RF_DATA/HYB_data/S_RF
else echo "Skiping"
fi
cd ./../
done
