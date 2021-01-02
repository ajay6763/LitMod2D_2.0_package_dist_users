#!/bin/bash

# THIS SCRIPT IS FOR DATA SAMPLED AT 20 SAMPLES PER SECOND.
# THAT IS MAINLY FOR IMD INDIAN NETWORK DATA "NT"
#
#
# station name where we have the events from.
# change this according to the station and network
# for which we are doing the receiver functions.
# 
# also change the delay accroding to the + and - cut from the t0 marker
# this is similar to the time offset in pwaveqn
#
# later make this code more general by putting this as a variable.

echo "Input 'revc' or 'revf' for unfiltered (c45) or filtered (f45) Data"
read pre

echo "Input the number of Iterations required (eg. 200)"
read nbumps
#nbumps=200

delay=45.0

# taking up all event directories starting with 200.......
for dir in [12]*
do

echo "working on $dir......."

# moving into the directory
cd $dir

# getting the station and network information from file name (using BHZ)
cmp=`ls ${pre}45_${dir}_*_*.?HZ | awk -F. '{print $2}' | cut -c1-2`
station=`ls ${pre}45_${dir}_*_*.?HZ | awk -F_ '{print $4}' | awk -F. '{print $1}'`
network=`ls ${pre}45_${dir}_*_*.?HZ | awk -F_ '{print $3}'`
echo "Station $station Network $network Component root $cmp"

# making directory for iterative deconvolution
mkdir Iterdecon_dir

# copying the cut 45 seconds file of the event, all components
cp ${pre}45_${dir}_${network}_${station}.${cmp}? Iterdecon_dir/

# moving into the Iterdecon_dir directory.....
cd Iterdecon_dir/

# decimating to 20 samples per second, 
# and rotating the horizontals to radial and transverse


sac << END
r ${pre}45_${dir}_${network}_${station}.${cmp}?
lh delta
r ${pre}45_${dir}_${network}_${station}.${cmp}N ${pre}45_${dir}_${network}_${station}.${cmp}E
rot to gcp
w ${pre}45_${dir}_${network}_${station}.${cmp}R ${pre}45_${dir}_${network}_${station}.${cmp}T
r ${pre}45_${dir}_${network}_${station}.${cmp}R
ch KCMPNM ${cmp}R
wh
r ${pre}45_${dir}_${network}_${station}.${cmp}T
ch KCMPNM ${cmp}T
wh
q
END


# listing all the c45_ files"
echo "files to be used are........."
ls ${pre}45_*.${cmp}Z ${pre}45_*.${cmp}R ${pre}45_*.${cmp}T
echo

# formulating the output file name. 

old_root=`echo $dir | cut -c3-11 `
echo "old root name is $old_root"
echo


# setting up the two different filters for calculating the receiver functions.

for water in 0.001 ;
do
echo " water level filter is $water"

for gauss in 1.5 2.5;
do
echo "Gaussian filter is $gauss"



# run a simple test receiver function example
#
#  note that the iterative vode is set up to
#    produce a result that when convolved with
#    the vertical, reproduces the radial. This
#    is not the case in the FDomain rftn. The
#    averaging function of the FDomain result
#    should have unit area, but it does not.
#    if you normalize the FDomain rftn, lac.eqr,
#    by the area in the rftn, you will get the
#    same amplitudes as the time domain.
#    See the macro compare.macro for details

# for the vertical component

iterdecon << end
${pre}45_${dir}_${network}_${station}.${cmp}Z
${pre}45_${dir}_${network}_${station}.${cmp}R
${nbumps}      * nbumps
${delay}      * phase delay for result
${water} * min error improvement to accept
${gauss}      * Gaussian width factor
1        * 1 allows negative bumps
0        * 0 form minimal output (1) will output lots of files
end
mv decon.out* ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqz

# Writing the component name in the SAC header
sac << END1
r ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqz
ch KCMPNM eqz
wh
q
END1


# for the radial component

iterdecon << end
${pre}45_${dir}_${network}_${station}.${cmp}R
${pre}45_${dir}_${network}_${station}.${cmp}R
${nbumps}      * nbumps
${delay}      * phase delay for result
${water} * min error improvement to accept
${gauss}      * Gaussian width factor
1        * 1 allows negative bumps
0        * 0 form minimal output (1) will output lots of files
end
mv decon.out* ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqr

# Writing the component name in the SAC header
sac << END2
r ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqr
ch KCMPNM eqr
wh
q
END2


# for the transverse component

iterdecon << end
${pre}45_${dir}_${network}_${station}.${cmp}T
${pre}45_${dir}_${network}_${station}.${cmp}R
${nbumps}      * nbumps
${delay}      * phase delay for result
${water} * min error improvement to accept
${gauss}      * Gaussian width factor
1        * 1 allows negative bumps
0        * 0 form minimal output (1) will output lots of files
end
mv decon.out* ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqt

# Writing the component name in the SAC header
sac << END3
r ${pre}45_${old_root}${station}_${water}_${gauss}.i.eqt
ch KCMPNM eqt
wh
q
END3


done

done

#rm *.${cmp}[RT] 
rm denominator* numerator* predicted* observed* name

cd ../..
echo "out of directory....moving to the next one"

done

echo ".......finished"
