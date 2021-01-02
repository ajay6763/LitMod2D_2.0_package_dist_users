#!/bin/bash

# Script to extract Receiver Function event list for a given station from NEIC-CSV list
# Script name: 

# Getting screen input of station name and location
echo "Station name"
read sta
echo "Station network"
read net

# Getting screen input of Catalogue file to search on
echo "Input the CSV Catalogue file to search (provide full path)"
echo "eg. /DATABACKUP/RESOURCES/CATLOGUES_DATABASES/NEIC_database/neic-csv_2013"
read catfile
# Else set it here
#catfile=/DATABACKUP/RESOURCES/CATLOGUES_DATABASES/NEIC_database/neic-csv_2013

# Getting the Start and End dates
#echo "Input year"
#read syear
#echo  "Input starting month"
#read smonth
#echo  "Input ending month"
#read emonth
# Else set it here
#syear=
#smonth=
#sday=
#emonth
#eday

# SORTING (RF events for the given STATION) ################
# Getting the Station lat long from STATION FILE
sta_file=/DATABACKUP/RESOURCES/Responses_and_Stations/ALL_station_Lists/all_station_file
stalat=`cat ${sta_file} | grep $net | grep $sta | awk '{print $4}'`
stalong=`cat ${sta_file} | grep $net | grep $sta | awk '{print $5}'`

echo "Extracting event list for Station $sta Lat $stalat Long $stalong"
# SORTING input catalogue file to zoom into the selected months and also magnitude 
# Magnitude 5.0 and above
#echo "Sorting by year, month, magnitude 5.0 ...."
cat $catfile |awk  '$11 >= 5 {print $0}'> temp5

rm -f cvsfile
cat temp5 > fish1

while [  -s fish1 ]; do
        line=`head -1 fish1`
        tail -n +2 fish1 > fish2
        cat fish2 > fish1
        rm -f fish2

#echo "line is $line"
evlat=`echo $line | awk  '{print $8}'`
evlong=`echo $line | awk  '{print $9}'`
# Calculating distance
dist=`/home/ajay/RF/Catlogue2eventlist/Catlogue2RF_eventlist/distaz_once << EOF | awk '{printf("%5.0f", $2)}' | awk '{print $1}'
$stalat $stalong
$evlat $evlong
EOF`
#echo "Distance is $dist"
if [ $dist -ge 30 ] && [ $dist -le 40 ]; then
echo "$line" >> cvsfile
fi

done

rm -f fish1 temp5 nul.tmp

cat $catfile |awk  '$11 >= 5.5 {print $0}'> temp5.5

cat temp5.5 > fish1

while [  -s fish1 ]; do
        line=`head -1 fish1`
        tail -n +2 fish1 > fish2
        cat fish2 > fish1
        rm -f fish2

#echo "line is $line"
evlat=`echo $line | awk  '{print $8}'`
evlong=`echo $line | awk  '{print $9}'`
# Calculating distance
dist=`/home/ajay/RF/Catlogue2eventlist/Catlogue2RF_eventlist/distaz_once << EOF | awk '{printf("%5.0f", $2)}' | awk '{print $1}'
$stalat $stalong
$evlat $evlong
EOF`
#echo "Distance is $dist"
if [ $dist -gt 40 ] && [ $dist -le 50 ]; then
echo "$line" >> cvsfile
fi

done

rm -f fish1 temp5.5 nul.tmp

cat $catfile |awk  '$11 >= 5.7 {print $0}'> temp5.7

cat temp5.7 > fish1

while [  -s fish1 ]; do
        line=`head -1 fish1`
        tail -n +2 fish1 > fish2
        cat fish2 > fish1
        rm -f fish2

#echo "line is $line"
evlat=`echo $line | awk  '{print $8}'`
evlong=`echo $line | awk  '{print $9}'`
# Calculating distance
dist=`/home/ajay/RF/Catlogue2eventlist/Catlogue2RF_eventlist/distaz_once << EOF | awk '{printf("%5.0f", $2)}' | awk '{print $1}'
$stalat $stalong
$evlat $evlong
EOF`
#echo "Distance is $dist"
if [ $dist -gt 50 ] && [ $dist -le 90 ]; then
echo "$line" >> cvsfile
fi

done

rm -f fish1 temp5.7 nul.tmp

mv cvsfile PALK.txt

# FORMATTING (selected events) ##############3

#outfile=${sta}_event.txt
# Writing the top line
#3echo "          Origin time           Lat      Lon     Dp  Mag" > ${outfile}

# Running through the input file
#cat cvsfile > fish1

#while [  -s fish1 ]; do

#        line=`head -1 fish1`
#        tail -n +2 fish1 > fish2
#        cat fish2 > fish1
#        rm -f fish2

#echo "line is $line"

#year=`echo $line | awk -F- '{print $1}'`
#month=`echo $line | awk -F- '{print $2}'`
#date=`echo $line | awk -F- '{print $3}' | awk -FT '{print $1}'`

#jdy=`~/bin/triton-bin/date_jdy << ! | tail -1
#$year
#$month
#$date
#!
#`
#hr=`echo $line | awk -F: '{print $1}' | awk -FT '{print $2}'`
#min=`echo $line | awk -F: '{print $2}'`
#sec=`echo $line | awk -F: '{print $3}' | awk -F. '{print $1}'`
#ss=`echo $line | awk -F: '{print $3}' | awk -F. '{print $2}' | awk -FZ '{print $1}'`#

#lat=`echo $line | awk -F, '{print $2}'`
#long=`echo $line | awk -F, '{print $3}'`
#depth=`echo $line | awk -F, '{print $4}'`
#mag=`echo $line | awk -F, '{print $5}'`

#id=${year}${jdy}${hr}${min}${sec}

# Writing the output in format
#echo "${year} (${jdy}) ${month} ${date} ${hr} ${min} ${sec} ${lat} ${long} ${depth} ${mag} ${id}" | awk '{printf("%4s %5s %2s %2s %2s %2s %4.1f %8.3f %8.3f %5.1f %4.1f %13s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12)}' >> ${outfile}

#done

# Cleaning up
rm -f fish1

