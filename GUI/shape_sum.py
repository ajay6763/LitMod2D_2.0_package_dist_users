from shapely.geometry.polygon import  Polygon as Shape
from shapely.ops import cascaded_union # to merge polygons
def shape_sum(x_envelope_1,x_envelope_2,y_envelope_1,y_envelope_2):
    a=Shape()
    poly_1=Shape()
    poly_2=Shape()
    poly_1=Shape( [(x_envelope_1[j],y_envelope_1[j])
    	for j in range(len(x_envelope_1)) ])
    poly_2=Shape( [(x_envelope_2[j],y_envelope_2[j]) 
    	for j in range(len(x_envelope_2)) ])
    
    polygons= [poly_1,poly_2]
    #a=poly_2.union(poly_1)
    a=cascaded_union(polygons)
    tx=[];ty=[];
    t_x,t_y=a.boundary.xy
    return t_x,t_y