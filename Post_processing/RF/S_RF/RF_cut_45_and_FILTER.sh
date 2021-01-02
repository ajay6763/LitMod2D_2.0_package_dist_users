#!/bin/bash

#	script name: RF_cut_45_and_FILTER.sh
#
#	This script 
#	 1. cuts the file 45 secs before and after the first p arrival 
#	 2. filters between 1 and 20 sec (on choice)
#	 3. writes it to a f45_ (filtered) or c45_ (unfiltered) SAC file
# 		These filtered/unfiltered RF data are used for profile plotting.
#

echo "Do you want to pre-filter data (Y/N)?"
read ans
if [ $ans = "y" ] || [ $ans = "Y" ] ; then
echo "Filtering and Processing..."
echo
# Filtering and chopping for RF
for dir in 20*
do
echo "event directory - $dir"
cd $dir
echo "into directory $dir....."
echo
	echo echo on > sac.mac
	echo r ${dir}*.?H? >> sac.mac
	echo qdp off >> sac.mac
	echo cut t0 -120 +45 >> sac.mac
	echo r >> sac.mac
	echo rmean >> sac.mac
	echo rtr >> sac.mac
	echo taper w 0.1 >> sac.mac	
	echo w prepend f45_ >> sac.mac
	echo r f45_* >> sac.mac
	echo hp co 0.05 p 2 n 2 >> sac.mac
	echo w over >> sac.mac
	echo ppk >> sac.mac
	echo quit >> sac.mac
${SACDIR}/bin/sac sac.mac
echo
rm -f sac.mac
# Going out of the directory and taking up the next one
cd ..
done
echo
echo "**** Remember to Decimate f45_ files to 20sps before Iterdecon *****"
echo

elif [ $ans = "n" ] || [ $ans = "N" ] ; then
echo "Processing without filtering..."
echo
for dir in 20*
do
echo "event directory - $dir"
cd $dir
echo "into directory $dir....."
echo
	echo echo on > sac.mac
	echo r ${dir}*.?H? >> sac.mac
	echo qdp off >> sac.mac
	echo cut t0 -120 +45 >> sac.mac
	echo r >> sac.mac
	echo rmean >> sac.mac
	echo taper w 0.1 >> sac.mac	
	echo w prepend c45_ >> sac.mac
	echo r c45_* >> sac.mac
	echo ppk >> sac.mac
	echo quit >> sac.mac
${SACDIR}/bin/sac sac.mac
echo
rm -f sac.mac
# Going out of the directory and taking up the next one
cd ..
done
echo
echo "**** Remember to Decimate c45_ files to 20sps before Iterdecon *****"
echo

else
echo "Answer mismatch. Should be Y or N. Quitting..."
exit
fi

# The end!
