#!/bin/bash
for dir in [12]*
do
cd $dir/Iterdecon_dir/
sachead 'f45_'$dir'_G_HYB.BHZ' user0 >ray_par
ray_para=`cat ray_par`
sac<<!
r *.eq?
ch user4 ${ray_para}
wh
q
!
rm ray_par
cd ../../

echo "moving to next directory"
done
