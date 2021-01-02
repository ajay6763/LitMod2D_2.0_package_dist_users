#!/bin/bash

# Initialising/cleaning up
rm -f logfile.txt percent_fit.txt

if [ -e best_eqrs ]; then 
echo "best_eqrs/ directory exisit ... cleaning up inside"
rm -f best_eqrs/*
else 
echo "Creating best_eqrs/"
mkdir best_eqrs
fi

# Best selection mode determination
echo "Do you want to select best eqrs by percent fit (f) or manually (m)"
read resp

# Automatic from percent fit
if [ $resp = "F" ] || [ $resp = "f" ] ; then

echo "Input the acceptable fit value as in integer (eg. 80, 85, 90 etc)"
read value

# Checking if the value entered is an integer
if echo $value | egrep -q '^[0-9]+$'; then 
echo "Value accepted. Proceeding"
else 
echo "Value should have been an integer (eg. 80, 85, 90 etc). Quitting..." 
exit
fi

# Running through all eqr files
for file in *eqr
do
# extracting the USER9 information from the sac header and rounding the value off to integer
fit=`saclhdr -USER9 $file | awk '{printf("%2.0f\n", $1)}'`
echo "$file  $fit" >> percent_fit.txt

if [ $fit -ge $value ]; then
echo "$file Fit = $fit selecting"
cp $file  best_eqrs/ 
elif [ $fit -lt $value ]; then
echo "$file Fit = $fit Not copying"
else
echo "Value mismatch. Quitting..."
exit
fi

done

# Manually selecting by viewing individual files
elif [ $resp = "M" ] || [ $resp = "m" ] ; then

for file in *.eqr
do

echo "working on eqr file: $file"

# Plotting data
${SACDIR}/bin/sac <<! 
r $file
qdp off
fileid on type list user9
xlim -40 120
ppk
q
!

echo
echo
echo "input 'y/n' for copying data or do nothing"
read ans

if [ ${ans} = "y" ]; then
echo "copying data to best directory"
cp $file  best_eqrs/
elif [ ${ans} = "n" ]; then
echo "doing nothing"
else 
echo "Response not understood: listing as not recognised"
echo "$file" >> logfile.txt
fi

# going on to the next file
done

else
echo "Answer not recognised. Should have been "f" or "m". Quitting..."
exit
fi

echo "Run this line in the best_eqrs directory to create a runme for copying eqts"
echo "ls f45_* | awk -F. '{print "cp /DATABACKUP/STUDENT_DATA/ajay_data/RF_DATA/HYB_data/HYB_IRIS_data/RF_eqt_2.5/"$1"."$2"."$3"."$4".eqt ."}' > runme"
echo "run the runme from the best_eqts directory"
# The end!
