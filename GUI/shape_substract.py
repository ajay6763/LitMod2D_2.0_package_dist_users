from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.ticker as ticker
from numpy import sqrt
import numpy as np
from shapely.geometry.polygon import  Polygon as Shape
from shapely.ops import cascaded_union # to merge polygons
from shapely.geometry import  Point, LinearRing
def shape_substract(x_envelope_1,x_envelope_2,y_envelope_1,y_envelope_2):
    a=Shape()
    poly_1=Shape()
    poly_2=Shape()
    poly_1=Shape( [(x_envelope_1[j],y_envelope_1[j])
    	for j in range(len(x_envelope_1)) ])
    poly_2=Shape( [(x_envelope_2[j],y_envelope_2[j]) 
    	for j in range(len(x_envelope_2)) ])
    a=poly_2.difference(poly_1)
    tx=[];ty=[];
    t_x,t_y=a.exterior.xy
    return t_x,t_y