#!/bin/bash

# This script calculates the piercing point of RFs using taup at a given depth
# Needs a control file to point to which station RFs to process called: station_file (see example_REPP/) 
# Runs through the RFs of the stations listed in the station_file

#echo "Input the RF list to use for computing the piercing points (eg. piercing_pt_input.txt)"
#read rflist

echo "Input the name of the output file to write to (eg. piercing_pt_taup.txt)"
read outfile

# Cleaning up
rm -f ${outfile}

for sta in `awk 'NR>1 {print $2}' station_file` ; do
moho=`grep $sta station_file | awk '{print $6}'`

# Running through all eqr files
for file in *${sta}*eqr ; do
# Getting header info for input to taup
stla=`saclhdr -STLA $file`
stlo=`saclhdr -STLO $file`
evla=`saclhdr -EVLA $file`
evlo=`saclhdr -EVLO $file`
evdp=`saclhdr -EVDP $file`

# Running TAUP (piercing point)
taup_pierce -h $evdp -ph P -sta  $stla $stlo  -evt $evla $evlo -o output.dat -nodiscon -pierce $moho

# checking output file length and deciding on the line to use
num=`wc -l output.dat | awk '{print $1}'`
if [ $num -eq 3 ] ; then 
plat=`awk 'NR==3{print $4}' < output.dat`
plong=`awk 'NR==3{print $5}' < output.dat`
else
plat=`awk 'NR==2{print $4}' < output.dat`
plong=`awk 'NR==2{print $5}' < output.dat`
fi
rm -f output.dat

# Writing piercing lat long to the output file
echo "$file $plong $plat" >> ${outfile}
# Looping over RF files
done
# Looping over stations
done

