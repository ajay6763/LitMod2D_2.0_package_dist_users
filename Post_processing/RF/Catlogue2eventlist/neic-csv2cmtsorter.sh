#!/bin/bash

echo "Input the name of the NEIC CSV formatted output file (from NEIC website)"
read csvfile

outfile=pde_event.txt
# Writing the top line
echo "          Origin time           Lat      Lon     Dp  Mag" > ${outfile}

# Running through the input file
tail -n +2 ${csvfile} > fish1

while [  -s fish1 ]; do

        line=`head -1 fish1`
        tail -n +2 fish1 > fish2
        cat fish2 > fish1
        rm -f fish2

echo "line is $line"

year=`echo $line | awk -F- '{print $1}'`
month=`echo $line | awk -F- '{print $2}'`
date=`echo $line | awk -F- '{print $3}' | awk -FT '{print $1}'`

jdy=`~/bin/triton-bin/date_jdy << ! | tail -1
$year
$month
$date
!
`
hr=`echo $line | awk -F: '{print $1}' | awk -FT '{print $2}'`
min=`echo $line | awk -F: '{print $2}'`
sec=`echo $line | awk -F: '{print $3}' | awk -F. '{print $1}'`
ss=`echo $line | awk -F: '{print $3}' | awk -F. '{print $2}' | awk -FZ '{print $1}'`

lat=`echo $line | awk -F, '{print $2}'`
long=`echo $line | awk -F, '{print $3}'`
depth=`echo $line | awk -F, '{print $4}'`
mag=`echo $line | awk -F, '{print $5}'`

id=${year}${jdy}${hr}${min}${ss}

# Writing the output in format
echo "${year} (${jdy}) ${month} ${date} ${hr} ${min} ${sec} ${lat} ${long} ${depth} ${mag} ${id}" | awk '{printf("%4s %5s %2s %2s %2s %2s %4.1f %8.3f %8.3f %5.1f %4.1f %13s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12)}' >> ${outfile}

done

# Cleaning up
rm -f fish1

