
rm -f all t data_phase_vel list
awk '{print $2}' 4s.dat   > data_phase_vel
for file in 6s.dat 8s.dat 10s.dat 12s.dat 14s.dat 16s.dat 18s.dat 20s.dat 25s.dat 30s.dat 35s.dat 40s.dat 45s.dat 50s.dat 59s.dat 67s.dat 77s.dat 87s.dat 100s.dat 111s.dat 125s.dat 143s.dat 167s.dat 
do 
echo $file
awk '{print $2}' $file > all
echo $file>> list
paste  data_phase_vel  all >> t 
mv t data_phase_vel
done



########### make dispersion curves along the profile
count=0
for i in `seq 0 5 670`
do
let "count=count+1"
cat data_phase_vel| awk -v var="$count" 'NR==var {print $0}'> disp_$i.dat
tr '\t' '\n' < disp_$i.dat > temp
mv temp disp_$i.dat
done
