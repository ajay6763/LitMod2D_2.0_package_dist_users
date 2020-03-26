from shapely.geometry.polygon import  Polygon as Shape
from shapely.ops import cascaded_union # to merge polygons
import shape_sum as shape_sum
def make(l,b):
    x_envelope=[]
    y_envelope=[]
    for i in range(len(l)):
        poly=[]
        poly_sum=[]
        x_envelope.append([])
        y_envelope.append([])
        tx=[]
        ty=[]
        for k in range(i,len(l)):
            tx,ty=shape_sum.shape_sum(l[k][:],l[k+1][:],b[k][:],b[k+1][:])
            poly.append(Shape([(l[k][j],b[k][j]) for j in range(len(l[k])) ]))
        poly_sum=cascaded_union(poly)
        x_envelope[i],y_envelope[i]=poly_sum.exterior.coords.xy
    x_envelope.append(l[len(l)-1])
    y_envelope.append(b[len(l)-1])
    return x_envelope,y_envelope