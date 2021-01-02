#!/bin/bash

mkdir problem_dir

for dir in 19* 20*
do 

num=`ls $dir | wc | awk '{print $1}'`

if [ $num -gt 3 ]; then 
echo "$dir is fine..."

else echo "$dir has data"
mv $dir problem_dir/

fi

done 
