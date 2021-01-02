#!/bin/bash

## Converts seisan event format to cmt format

# cleaning up
rm -f naman_event.txt tmp1 tmp2

# Running through the list
echo "Input the list of seisan files to be used (eg. list: created by ls 18* > list)"
read filelist

cat ${filelist} > fish1

while [  -s fish1 ]; do

        line=`head -1 fish1`
        tail -n +2 fish1 > fish2
        cat fish2 > fish1
        rm -f fish2

file=`echo $line | awk '{print $1}'`
echo
echo "Working on file $file"

# Checking the month as 10 and above is written with a space in the input file
monthcheck=`ls $file | awk -F. '{print $2}' | cut -c6-7`
echo $monthcheck
if [ $monthcheck -ge 10 ]; 
then
year=`awk 'NR==1 {print $1}' $file` 
monthi=`awk 'NR==1 {print $2}' $file`
month=`awk -F- 'NR==3 {print $2}' $file`
datei=`awk 'NR==1 {print $3}' $file`
date=`awk -F- 'NR==3 {print $3}' $file`

else
year=`awk 'NR==1 {print $1}' $file` 
monthi=`awk 'NR==1 {print $2}' $file | cut -c1`
month=`awk -F- 'NR==3 {print $2}' $file`
datei=`awk 'NR==1 {print $3}' $file | cut -c2-3`
date=`awk -F- 'NR==3 {print $3}' $file`
fi

hr=`awk -F- 'NR==3 {print $4}' $file | cut -c1-2`
min=`awk -F- 'NR==3 {print $4}' $file | cut -c 3-4`
ss=`awk 'NR==1 {print $5}' $file | awk -F. '{print $2}'`
sec=`awk -F- 'NR==3 {print $5}' $file | awk -FS '{print $1}'`

echo "Year $year month $month date $date hour $hr minute $min second ${sec}.${ss}"

# Calculate jdy
/home/mitra/bin/date2jdy.py <<! > tmp1
$year
$monthi
$datei
!

jdy=`awk 'NR==2 {print $4}' tmp1`

	if [ $jdy -lt 10 ]; then
        jdy="00"$jdy
         elif [ $jdy -gt 10 ] && [ $jdy -lt 100 ]; then
         jdy="0"$jdy
         fi
#jdy=000

echo "$year ($jdy) $month $date $hr $min ${sec}.${ss} ${year}${jdy}${hr}${min}${sec}"
echo "$year ($jdy) $month $date $hr $min ${sec}.${ss} ${year}${jdy}${hr}${min}${sec}" >> tmp2

done

# Writing the event file in cmt format
echo "          Origin time           Lat      Lon     Dp  Mag" >  naman_event.txt
awk '{printf("%4s %5s %2s %2s %2s %2s %4.1f %8.3f %9.3f %4s %4.1f %13s\n", $1, $2, $3, $4, $5, $6, $7, "00.000", "00.000", "00", "0.0", $8)}' tmp2 >> naman_event.txt

# cleaning up
rm -f fish1 tmp1 tmp2
