"""
This is an example to show how to build cross-GUI applications using
matplotlib event handling to interact with objects on the canvas
####
"""
from matplotlib.artist import Artist
#from matplotlib import Polygon, CirclePolygon
from numpy import sqrt, nonzero, equal, array, asarray, dot, amin, cos, sin,searchsorted
import numpy as np
from matplotlib.mlab import dist_point_to_segment

from matplotlib.patches import Polygon
from matplotlib.pyplot import *
import matplotlib as plt

class PolygonInteractor:
    """
    An polygon editor.

    Key-bindings

      't' toggle vertex markers on and off.  When vertex markers are on,
          you can move them, delete them

      'd' delete the vertex under point

      'i' insert a vertex at point.  You must be within epsilon of the
          line connecting two existing vertices

    """

    showverts = True
    epsilon = 5  # max  distance to count as a vertex hit

    def __init__(self, ax_edit, poly,no_nodes,digitized_file,reso,Z,topo,out):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure or canvas before defining the interactor')
        self.ax_edit = ax_edit
        self.out_file=out;
        canvas = poly.figure.canvas
        self.topo=topo
        self.poly = poly
        self.l=[]
        self.b=[]
        self.reso=reso
        self.Z=Z
        self.cum_l=[]
        self.no_nodes=no_nodes
        self.cum_l=np.ndarray.tolist(np.cumsum(self.no_nodes))
        self.cum_l.insert(0,0)
        print('Index counter to retreive geometries')
        print self.cum_l
        x, y = zip(*self.poly.xy)
        self.line = Line2D(x,y,marker='o',markersize=4, markerfacecolor='r', animated=True)
        self.ax_edit.add_line(self.line)
        cid = self.poly.add_callback(self.poly_changed)
        print('length of the poly')
        #print len(self.poly)
        print('Sum of the number of nodes')
        print sum(self.no_nodes)
        self._ind = None # the active vert
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_press_event', self.save_poly)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax_edit.bbox)
        #self.ax.draw_artist(self.poly)
        self.ax_edit.draw_artist(self.line)
        try:
            f=np.loadtxt(digitized_file,comments=">")
            self.ax_edit.plot(f[:,0],f[:,1],'*',lw=1,color='grey')
        except:
            pass
        try:
            f=np.loadtxt(self.topo)
            self.ax_edit.plot(f[:,0],f[:,1]/1000.0,color="black",lw=0.1,marker='d',markersize=5)
        except:
            pass
        self.canvas.blit(self.ax_edit.bbox)

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, self.poly)
        self.line.set_visible(vis)  # don't use the poly visibility state
    def save_poly(self,event):
    	if event.button == 3:
    		self.l=self.poly.xy[:,0]
    		self.b=self.poly.xy[:,1]
            ### fihure the cumulative lengths of the envelops
    		#print "printing modofied coordinates"
    		#print self.l
    		#print self.b
    		print "total length"
    		print len(self.l)
    		f=open(self.out_file,"w")
    		f.writelines(">  \t  %d\n " %(len(self.no_nodes)))
    		# this will take care of i=0 problem
    		#f.writelines(">  \t  %d\n " % ((self.no_nodes[0])))
    		j=-1
    		for i in range(len(self.no_nodes)):
    			#if i==0:
    			#	j=j+1
    			#	f.writelines(">  \t  %d\n " % ((self.no_nodes[j])))
    		#	print i,self.no_nodes[j]
    		#	if i==self.no_nodes[j]:
    		#		print "new body"
    		#		j=j+1
    		#		f.writelines(">  \t  %d\n " % ((self.no_nodes[j])))
    		#		f.writelines(["%f \t %f \n" % (self.l[i] , self.b[i])])
    		#	else:
    			f.writelines(">  \t  %d\n " % ((self.no_nodes[i])))
    			for j in range(self.cum_l[i],self.cum_l[i+1]):
    				f.writelines(["%f \t %f \n" % (self.l[j],self.b[j])])
    		f.close() 
    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        # display coords
        xy = asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = nonzero(equal(d, amin(d)))[0]
        ind = indseq
        print ind
        if d[ind[0]]>=self.epsilon:
            ind = None
        return ind
    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)
    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None
    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        # part where insertion and deletion od points is done. Commenting this for now will come back to this later
        if event.key=='t':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts: self._ind = None
        elif event.key=='d':
            ind = self.get_ind_under_point(event)
            self.index=ind
            #if len(ind)==1:
            if ind is not None:
                for i in range(len(self.index)):
                    for kk in range(len(self.index)):
                        ind = self.get_ind_under_point(event)
                        self.index=ind
                        self.poly.xy = [tup for i,tup in enumerate(self.poly.xy) if i!=self.index[0]]
                        self.line.set_data(zip(*self.poly.xy))
                        ## figuring in which body insertion is done
                        print('Nodes before deletion')
                        print self.cum_l
                        here=searchsorted(self.cum_l,self.index[0])
                        print here
                        print('Index')
                        print self.index[0]
                        for jj in range(here,len(self.cum_l)):
                            print jj
                            self.cum_l[jj]=self.cum_l[jj]-1
                        print('Nodes after deletion')
                        print self.cum_l
                        break
        elif event.key=='i':
            print('Index at which insertion is done')
            ind = self.get_ind_under_point(event)
            print('No of nodes before insertion')
            print(self.cum_l)
            self.index=ind
            if self.index.all() != None:
                for kk in range(len(self.index)):
                    ind = self.get_ind_under_point(event)
                    self.index=ind
                    #print('Indices where to insert')
                    #print self.index[:]
                    self.x=np.ndarray.tolist(self.poly.xy[:,0])
                    self.y=np.ndarray.tolist(self.poly.xy[:,1])
                    if self.index[kk]>0:
                        #print('insertion at index')
                        #print self.index[kk]
                        #print('inserted values')
                        self.poly.xy = array(
                            list(self.poly.xy[:self.index[kk]]) +
                            [(self.x[self.index[0]]-self.reso, self.y[self.index[0]])] +
                            list(self.poly.xy[self.index[kk]:]))
                        self.line.set_data(zip(*self.poly.xy))
                        here=searchsorted(self.cum_l,self.index[kk])
                        print('Here will be added')
                        print self.cum_l[here]
                        #print('Index')
                        #print self.index[00]
                        for jj in range(here,len(self.cum_l)):
                            print jj
                            self.cum_l[jj]=self.cum_l[jj]+1
                        print('No of nodes after' + str(kk) + ' insertion =' + str(self.cum_l))
                            #print(self.cum_l)
                self.index=None
            else:
                pass 
        self.canvas.draw()
        return
    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        #x =round(event.xdata/float(self.reso))*float(self.reso)
        #x =round(event.xdata/1.0)*1.0
        #y=round(event.ydata)
        #y = event.ydata
        x,y = event.xdata, event.ydata
        x=round(x)
        self.poly.xy[self._ind] = x,y
        #print self.poly.xy
        #self.l, self.b = zip(*self.poly.xy)
        self.line.set_data(zip(*self.poly.xy))
        self.canvas.restore_region(self.background)
        #self.ax.draw_artist(self.poly)
        self.ax_edit.draw_artist(self.line)
        self.canvas.blit(self.ax_edit.bbox)


