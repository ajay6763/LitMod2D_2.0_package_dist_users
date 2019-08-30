#!/bin/sh
P=0.06
for file in *.dat
do
for ALP in 2.5
do
hrftn96 -P -ALP ${ALP} -DT 0.1 -D 5. -RAYP ${P} -M $file  -2 hr	-NSAMP 1500	
mv hrftn96.sac $file.eqr
done
done
