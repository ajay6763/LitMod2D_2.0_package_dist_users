#!/bin/sh
P=0.08
for CALP in 25
do
case ${CALP} in
25) ALP=2.5 ;;
esac
hrftn96 -P -ALP ${ALP} -DT 0.1 -D 10. -RAYP ${P} -M model.dat
mv hrftn96.sac ${CALP}.rfn
done
