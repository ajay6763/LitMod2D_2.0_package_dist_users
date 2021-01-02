#!/bin/bash
for dir in [12]*
do
cd $dir/Iterdecon_dir
cp f45*.001_0.6.i.eq[rt] ../../RF_0.6
cd ../../
done
