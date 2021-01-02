#!/bin/bash
#####################################################################
#################	September 2019	@ajay	#####################
#####################################################################
##### This script works for upto 4th version of GMT             #####
##### Task- need to adapt this for GMT5 or higher version       #####
#####If you do show please send me a copy at ajay6763@gmail.com #####
#####################################################################
#####   Extract dispersion curves from global GDM52 Ekström (2011)###
#####Input: 							#####
#####		1. start and end Lat,Long of the profiled	#####
#####		2. Resolution used in LitMod (2km,5km..)	#####
#####Output:	Dispersion curves at each distant point along	##### 
#####           the profile  					#####
#####################################################################
##Note:
## The fourth column is the deviation with respect to PREM

#####################################
#### gettin some info about model
pwd 
echo 'Enter the starting Long,Lat of your profile:'
#read A
## Alboran
#A="-4.3/38.9"
## Algerian
A="1.0/41.8"
echo
echo 'Enter the ending Long,Lat of your profile:'
#read B
## Alboran
#B="-1.985/33.59"
## Algerian
B="5.2375/34.8635"
echo
echo 'Enter the resolution used in your profile (e.g. 2km,5km.. etc):'
read reso

###########################################
#### Now will get the Lat, Long along the profile with resolution
gridfile=/home/akumar/owncloud/PHD/RESOURCES/GMT/GRD_FILES/mediterranean.grd
project  -C$A -E$B -G$reso -Q > Long_Lat_along_profile.xy



####################################
### Now, each point along the profile will be read
### Disperion curves will be extracted and will be saved 
### with distance in Km in the name of the output file

ans=1
echo
echo 'What kind of sufrace wave you want (Rayleigh---1 (default case) , Love--0):'
read ans
surf_type=R
if [ $ans -eq 0 ]
then
surf_type=L
while read line; do
Long=`echo "$line"|awk '{print $1}'`
Lat=`echo "$line"|awk '{print $2}'`
dist=`echo "$line"|awk '{print $3}'`
#############################
### for Rayleigh waves
./GDM52_dispersion <<!
2
$Lat $Long
99
!
#mv GDM52_dispersion.out $dist"km_L".txt
########################################################
### Making input file for surf96
### Love wave phase velocities
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","L  C  X  0", 1000/($1), $3,($4*$3)/100)}'> temp
tac temp > $dist"km_L_C_SURF96".inp
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","L  C  X  0", 1000/($1), $2,0)}'> temp 
tac temp > PREM_L_C_SURF96.inp
### Love wave group velocities
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","L  U  X  0", 1000/($1), $6,($6*$7)/100)}'> temp #$dist"km_L_U_SURF96".inp
tac temp> $dist"km_L_U_SURF96".inp
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","L  U  X  0", 1000/($1), $5,0)}'> temp #PREM_L_U_SURF96.inp
tac temp > PREM_L_U_SURF96.inp
#rm  GDM52_dispersion.out 
done <Long_Lat_along_profile.xy
else
surf_type=L
while read line; do
Long=`echo "$line"|awk '{print $1}'`
Lat=`echo "$line"|awk '{print $2}'`
dist=`echo "$line"|awk '{print $3}'`
#############################
### for Rayleigh waves
./GDM52_dispersion <<!
77
2
$Lat $Long
99
!
#mv GDM52_dispersion.out $dist"km_R".txt
########################################################
### Making input file for surf96
### Rayleigh wave phase velocities
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","R  C  X  0", 1000/($1), $3,($4*$3)/100)}'> temp #$dist"km_R_C_SURF96".inp
tac temp > $dist"km_R_C_SURF96".inp
rm -f temp
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","R  C  X  0", 1000/($1), $2,0)}'> temp 
tac temp > PREM_R_C_SURF96.inp
rm -f temp
### Rayleigh wave group velocities
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","R  U  X  0", 1000/($1), $6,($6*$7)/100)}'> temp 
tac temp > $dist"km_R_U_SURF96".inp
rm -f temp
tail -n +6 GDM52_dispersion.out |awk '{printf("%s\t%s\t%f\t%f\t%f\n" ,"SURF96","R  U  X  0", 1000/($1), $5,0)}'> temp
tac temp >PREM_R_U_SURF96.inp
rm -f temp
#rm  GDM52_dispersion.out 
done <Long_Lat_along_profile.xy
fi
