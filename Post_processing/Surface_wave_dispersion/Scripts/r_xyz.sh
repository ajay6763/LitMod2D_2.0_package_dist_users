#/!bin/bash
for file in *.eqr
do
sac2xy $file $file.xyz
done
