#!/bin/bash

# script name: put_t0_timemark.sh
#
#				(NOT generalised set the year)
#
# Puts the t0 timemark at the P-peak in the BB data, for clarity of time 
# picking the data is bandpass filtered between 1 and 100 seconds period.

# taking up all directories starting with 19?? and 200?
# change this according to the year you are working on
#
for dir in 2*
do


# moving into directiry to plot BHZ.

cd $dir


# plotting the receiver functions for viewing
$SACDIR/bin/sac <<!
r *?H?
qdp off
rtr
rmean
bp co .01 1 p 2 n 3
ppk
wh
q
!


# getting the t0 timemark time
time=`$SACDIR/bin/sac <<! |grep "t0" | awk '{ print $3 }'
r *?HN
lh t0
q
!`

# marking t0 on all components
$SACDIR/bin/sac <<!
qdp off
r *?H?
ch t0 ${time}
wh
r
ppk
q
!


cd ..


# loop for the directory ends here, will move to the next directory.

done


# the end
