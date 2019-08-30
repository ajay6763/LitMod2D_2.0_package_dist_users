import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
 
data=np.loadtxt("disp_all.out")
data_=np.loadtxt("disp_all.inp")

colors = ("red", "green", "blue")
groups = ("coffee", "tea", "water") 
# Create plot
fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1, axisbg="1.0")
ax = fig.gca()

plt.imshow(data[:,0],data[:,1],data_[:,2]-data[:,2])
plt.colorbar()
plt.show()
'''
ax.scatter(data[:,0],data[:,1],data_[:,2]-data[:,2], alpha=0.8, edgecolors='none', s=30)
#ax.set_zlim(2.0,5.0)
plt.title('Matplot 3d scatter plot')
plt.legend(loc=2)
plt.show()
'''