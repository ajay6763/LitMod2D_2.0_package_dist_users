from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.ticker as ticker
from numpy import sqrt
import numpy as np
from shapely.geometry.polygon import  Polygon as Shape
from shapely.ops import cascaded_union # to merge polygons
from shapely.geometry import  Point, LinearRing
import math

def shape_split(vert,l_split,b_split,l,b):
    x_drawn = [vert[k][0] for k in range(len(vert))]
    y_drawn = [vert[k][1] for k in range(len(vert))]
    x_st=x_drawn[0] #x_st[n_bod-1]=x
    y_st=y_drawn[0]#y_st[n_bod-1]=y
    x_end=x_drawn[-1]
    y_end=y_drawn[-1]
    ## figuring out where does this fall in body to be split
    print ("start index for split")
    temp_d=[]
    temp_d=[ math.sqrt((l_split[i] - x_st)**2 + 
        (b_split[i] - y_st)**2) for i in range(len(l_split))]
    split_st_index=temp_d.index(min(temp_d))
    print ("end index for print ")
    temp_d=[]
    temp_d=[ math.sqrt((l_split[i] - x_end)**2 + 
        (b_split[i] - y_end)**2) for i in range(len(l_split))]
    split_end_index=temp_d.index(min(temp_d))
    print('updating path ')
    path=[]
    vert = []
    print ("adding in the start from the actual body")
    for k in range(0,split_st_index):
        print (k)
        vert.append((l_split[k],b_split[k]))
    print ("adding the drwan points")
    for f in range(len(x_drawn)):
        vert.append((x_drawn[f],y_drawn[f]))
    #    path, = ax.plot([x_drawn[f]],[y_drawn[f]],'o-',lw=1,color=color)
    # just append the end part of the split axis buy but have append what ever is drwan
    # append the end point of the axis
    print ("adding in the end from the actual body") 
    for k in range(split_end_index,len(l_split)):
        print (k)
        vert.append((l_split[k],b_split[k]))
    
    #    path, = ax.plot([l_split[k]],[b_split[k]],'o-',lw=1,color=color)
    # close the body
    ## calculating inserted body 
    temp_1_shape=Shape(vert[:][:])
    temp_2_shape=Shape([( l[k], b[k]) 
        for k in range(len(l))])
    a=Shape();xa=[];ya=[]
    a=temp_2_shape.difference(temp_1_shape)# this the left over part of the splitted body
    err=a.is_valid
    print (a.exterior.coords)
    xa,ya=a.exterior.xy                  
    return xa,ya,a,err
