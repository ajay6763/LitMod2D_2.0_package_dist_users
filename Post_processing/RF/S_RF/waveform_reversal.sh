#!/bin/bash
## script to reverse the waveform using a small c program written by me
## It needs a list of files to be convetred
echo "Input 'c' or 'f' for unfiltered (c45) or filtered (f45) Data"
read pre
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
ls ${pre}45_${dir}_${network}_${station}.${cmp}Z > t
sac2xy<t
rm t
/home/ajay/RF/S_RF/waveform_reversal << end
${pre}45_${dir}_${network}_${station}.${cmp}Z.xy 
end
xy2sac out.txt 'rev'${pre}45_${dir}_${network}_${station}.${cmp}'Z'
rm out.txt


ls ${pre}45_${dir}_${network}_${station}.${cmp}E > t
sac2xy<t
rm t
/home/ajay/RF/S_RF/waveform_reversal << end
${pre}45_${dir}_${network}_${station}.${cmp}E.xy 
end
xy2sac out.txt 'rev'${pre}45_${dir}_${network}_${station}.${cmp}'E'
rm out.txt

ls ${pre}45_${dir}_${network}_${station}.${cmp}N > t
sac2xy<t
rm t
/home/ajay/RF/S_RF/waveform_reversal << end
${pre}45_${dir}_${network}_${station}.${cmp}N.xy 
end
xy2sac out.txt 'rev'${pre}45_${dir}_${network}_${station}.${cmp}'N'
rm out.txt
## feeding header info
# station latitudesachead 'f45_'$dir'_G_HYB.BHZ' user0 >ray_par

file=${pre}45_${dir}_${network}_${station}.${cmp}Z;
STLA=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z STLA`;

EVLA=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z EVLA`;

STLO=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z STLO`;

STEL=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z STEL`;

EVLO=` sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z EVLO`;

EVDP=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z EVDP`;

DIST=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z DIST`;
AZ=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z AZ`;
BAZ=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z BAZ`;
GCARC=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z GCARC`;
KZDATE=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z KZDATE`;
KZTIME=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z KZTIME`;
KSTNM=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z KSTNM`;

KNEWK=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z KNETWK`

USER0=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z USER0`;

USER7=`sachead ${pre}45_${dir}_${network}_${station}.${cmp}Z USER7`;
echo echo on > sac.mac
echo r rev*  >> sac.mac
echo ch STLA ${STLA} >>sac.mac
echo ch EVLA ${EVLA} >>sac.mac
echo ch STLO ${STLO} >>sac.mac
echo ch EVLO ${EVLO} >>sac.mac
echo ch STEL ${STEL} >>sac.mac
echo ch EVDP ${EVDP} >>sac.mac
echo ch DIST ${DIST} >>sac.mac
echo ch AZ ${AZ} >>sac.mac
echo ch BAZ ${BAZ} >>sac.mac
echo ch GCARC ${GCARC} >>sac.mac
echo ch KSTNM $KSTNM >>sac.mac
#echo ch KZTIME $KZTIME >>sac.mac
echo ch KNETWK $KNEWK >>sac.mac
echo ch USER0 ${USER0} >>sac.mac
echo ch USER7 ${USER7} >>sac.mac
echo w over >>sac.mac
echo quit >>sac.mac
sac sac.mac
sac << !
r rev${pre}45_${dir}_${network}_${station}.${cmp}Z
ch CMPAZ 0
ch CMPINC 0
ch KCMPNM BHZ
wh
q
!
sac << !
r rev${pre}45_${dir}_${network}_${station}.${cmp}E
ch CMPAZ 90
ch CMPINC 90
ch KCMPNM BHE
wh
q
!
sac << !
r rev${pre}45_${dir}_${network}_${station}.${cmp}N
ch CMPAZ 0
ch CMPINC 90
ch KCMPNM BHN
wh
q
!
cd ./../
done
