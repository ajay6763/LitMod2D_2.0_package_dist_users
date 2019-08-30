#!/bin/bash



##### get lat long along the profile

## Alboran profile
A="-4.8/38.7" # "-4.0/39"   #"-4.3/38.9"
B="-2/32.8" #"-2/32.8"   #"-1.7/32.9"

## Algerean profile
#A="1.0/41.8"
#B="5.3/34.75"



slat=`echo $A | awk -F"/" '{print $2}'`
slong=`echo $A | awk -F"/" '{print $1}'`
elat=`echo $B | awk -F"/" '{print $2}'`
elong=`echo $B | awk -F"/" '{print $1}'`
topo=/home/akumar/owncloud/PHD/RESOURCES/GMT/GRD_FILES/mediterranean.grd
# Plotting the topography in the profile

# getting the totoal length of the profile from topo.xz file


project $topo -C$A -E$B -G5 -Q > t1.xy 


for file in *.*.eqr
do
i=`echo $file| awk -F_ '{print $1}'`
echo "line number"
echo $i
awk -v var=$i 'FNR==var  {print $0}'  t1.xy > temp
more temp

long=`awk '{print $1}' temp `
lat=`awk '{print $2}' temp `
# marking t0 on all components
sac <<!
r $file
ch STLO $long
ch STLA $lat
ch USER4  0.06
ch BAZ 90
wh
r
q
!

done
