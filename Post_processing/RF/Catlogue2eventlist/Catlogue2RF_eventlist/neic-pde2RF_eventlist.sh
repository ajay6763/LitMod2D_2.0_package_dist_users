#!/bin/bash
# Script to run CMT_sorter to obtain Receiver Function event list for a given station
# Script name: 

# Getting screen input of station name and location
echo "Station name"
read sta
echo "Input station laitude"
read stalat
echo "Input station longitude"
read stalong

# setting the start and end time and date
sdate=20120601
stime=0000
edate=20121231
etime=2400

# setting the catlogue file to search
#catfile=/DATABACKUP/CATLOGUES_DATABASES/NEIC_database/pde_2011
catfile=/DATABACKUP/RESOURCES/CATLOGUES_DATABASES/NEIC_database/pde_2012

# Running CMT sorter to extract RF event list from NEIC catlog on SERVER
# GCARC dist 30-40 degrees; Magnitude 5.0 and above 
CMT_sorter -cat pde $catfile -cir $stalat $stalong 30. 40. -mag 5.0 10.0 -date $sdate $stime $edate $etime -jdy
cat search.txt > ${sta}_event.txt

# GCARC dist 40-50 degrees; Magnitude 5.5 and above
CMT_sorter -cat pde $catfile -cir $stalat $stalong 40. 50. -mag 5.5 10.0 -date $sdate $stime $edate $etime -jdy
cat search.txt >> ${sta}_event.txt

# GCARC dist 50-90 degrees; Magnitude 5.7 and above
CMT_sorter -cat pde $catfile -cir $stalat $stalong 50. 90. -mag 5.7 10.0 -date $sdate $stime $edate $etime -jdy
cat search.txt >> ${sta}_event.txt

# Sorting by date
sort ${sta}_event.txt > temp
mv temp ${sta}_event.txt

