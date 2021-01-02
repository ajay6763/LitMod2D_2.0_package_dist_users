echo on
r $1
rmean
rtrend
ch o
ch T0
ch T1
ch T2
ch T3
ch T4
ch T5
ch T6
ch T7
ch T8
ch T9
w over
r $1
ch evla 39.900 evlo 142.020 evdp 15
evaluate to msec 5 * 100
ch o gmt 2012 087 11 00 44 %msec
w over
r $1
setbb otime &1,o
evaluate to tshift -1 * %otime
ch allt %tshift
ch user7 6.0
ch lovrok true
w over
cut o 0 3600
cuterr fillz
r
cut off
rmean
w over
quit
