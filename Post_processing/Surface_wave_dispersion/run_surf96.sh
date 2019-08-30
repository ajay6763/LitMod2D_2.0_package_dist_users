#!/bin/bash 
#cd ../cluster_1
##### 
# 	clean up 
##### 
surf96 39 

##### 
# 	define damping 
##### 
#surf96 32 0.1

##### 
# 	Select differential smoothing 
##### 
#surf96 36 1 

##### 
# 	Invert for Vp and Vs keeping Vp/Vs fixed 
##### 
#surf96 30 1 

#################################################
# changing weight on the 2 crustal layers
#surf96 45 
#surf96 31 10 0 

# checking
#surf96 45

##### 
#	set up repeated run for 5 iterations 
##### 
# 1 iteration: for forward modeling and tests
#surf96 1 2 6 1 2 

# 5 iteration:
#surf96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2

#20 Iterations 
surf96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 

#30 Iterations 
surf96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 

##### 
# plot the model and show the data fit after 5 iterations 
##### 
srfphv96 
plotnps -EPS -K -F7 -W10 < SRFPHV96.PLT > figsrf1.eps 

##### 
# 	write the output dispersion curve 
#####
surf96 27 disp.out

##### 
# 	save current model 
##### 
surf96 28 modl.out 

###########################################
#       plot resolution kernel

#plot resolution kernel on screen
#surf96 9

surf96 29 res-kernl.txt

srfphr96
mv SRFPHR96.PLT R.PLT
plotnps -EPS -K < R.PLT > figsrfn.eps
##########################################

##### 
# 	compare the individual models from the inversion 
# 	to the true model 
##### 
shwmod96 -K 1 -W 0.05 model.true 
mv SHWMOD96.PLT T.PLT 

shwmod96 -K -1 tmpmod96.??? 
mv SHWMOD96.PLT I.PLT 

cat T.PLT I.PLT > IT.PLT 
plotnps -EPS -K -F7 -W10 < IT.PLT > figsrf2.eps

gv figsrf1.eps

# The end
