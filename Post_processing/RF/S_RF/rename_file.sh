#!/bin/bash

echo "Input 'c' or 'f' for unfiltered (c45) or filtered (f45) Data"
read pre


# taking up all event directories starting with 200.......
for dir in 20*
do

echo "working on $dir......."

# moving into the directory
cd $dir

# getting the station and network information from file name (using BHZ)
cmp=`ls ${pre}45_${dir}_*_*.?HZ | awk -F. '{print $2}' | cut -c1-2`
station=`ls ${pre}45_${dir}_*_*.?HZ | awk -F_ '{print $4}' | awk -F. '{print $1}'`
network=`ls ${pre}45_${dir}_*_*.?HZ | awk -F_ '{print $3}'`
echo "Station $station Network $network Component root $cmp"

mv ${pre}45_${dir}_${network}_${station}_.${cmp}Z ${pre}45_${dir}_${network}_${station}.${cmp}Z
mv ${pre}45_${dir}_${network}_${station}_.${cmp}N ${pre}45_${dir}_${network}_${station}.${cmp}N
mv ${pre}45_${dir}_${network}_${station}_.${cmp}E ${pre}45_${dir}_${network}_${station}.${cmp}E

cd ..
done
