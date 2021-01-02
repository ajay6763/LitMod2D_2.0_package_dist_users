#!/bin/bash

rm -f logfile.txt
#mkdir best_eqrs_by_match_per

for file in *_0.6.i.eqr 
do
echo "working on eqr file: $file"
echo $file

per=`sachead $file USER9`
echo $per
l_b=90;
if [ 1 -eq `echo "$per >= $l_b" | bc` ]; then
echo "copying data to best directory"
cp $file  best_eqrs_by_match_per/
else
echo "Doing nothing"
fi
echo
# going on to the next file
done

echo "Run this line in the best_eqrs directory to create a runme for copying eqts"
echo "ls f45_* | awk -F. '{print "cp /DATABACKUP/STUDENT_DATA/ajay_data/RF_DATA/HYB_data/HYB_IRIS_data/RF_eqr_0.6/"$1"."$2"."$3"."$4".eqt ."}' > runme"
echo "run the runme from the best_eqts directory"
# The end!
