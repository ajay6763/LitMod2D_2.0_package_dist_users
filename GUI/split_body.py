from matplotlib import pyplot as plt
from matplotlib.widgets import Button, Cursor
from matplotlib.patches import Polygon


class Split_axis:
    def __init__(self,fig,ax,done_option):
        self.ax=ax
        self.fig=fig
        self.done_option=done_option
        self.canvas=self.fig.canvas
        self.tx =[] #list(line.get_xdata())
        self.ty =[] #list(line.get_ydata())
        #self.cid = line.figure.canvas.mpl_connect('button_press_event', self)
        self.canvas.mpl_connect('button_press_event', self.select_point)
        self.mouse_button = {1: self.select_point}

    def select_point(self, event):
        #print('click', event)
        #if event.inaxes!=self.line.axes: return
        if event.inaxes!=self.ax: return
        if event.button==1:
            #self.xs=event.xdata
            #self.xy=event.ydata
            self.tx.append(event.xdata)
            self.ty.append(event.ydata)
            print event.xdata,event.ydata
            #self.line.set_data(self.xs, self.ys)
            self.ax.plot(self.tx,self.ty)
            #self.line.figure.canvas.draw()
            plt.show()
    def done(self,event):
        if event.inaxes==self.done_option:
            print "quiting"
            f=open("split_body_dump.dat","w")
            for i in range(len(self.tx)):
                f.writelines(["%f \t %f \n" % (self.tx[i],self.ty[i])])
            plt.close(self.fig)

def insert_body(l_split,b_split,X,Y,major_xticks,major_yticks):
    fig_split = plt.figure(1)
    fig_split.canvas.set_window_title('Litmod Split Body')
    ax_split = fig_split.add_subplot(111)
    cursor = Cursor(ax_split, useblit=True, color='white', linewidth=0.5)
    plt.subplots_adjust(bottom=0.2)
    done_option = plt.axes([0.1, 0.05, 0.1, 0.075])
    ax_split.set_xlim((-2, 2))
    ax_split.set_ylim((-2, 2))
    ax_split.set_xlim((0, X))
    ax_split.set_ylim((-Y,0))
    ax_split.set_yticks(major_yticks)
    ax_split.set_xticks(major_xticks)  
    ax_split.grid(True)
    poly = Polygon(list(zip(l_split,b_split)),lw=1)
    ax_split.add_patch(poly)
    #ax_split.plot(l_split,b_split,marker='o', markerfacecolor='r', animated=True)
    p=Split_axis(fig_split,ax_split,done_option)
    done = Button(done_option, 'done',color="black")
    done.on_clicked(p.done)
    
    print "axis along which to split"
    print p.tx,p.ty
    plt.show()
    exit()
    


"""
      
fig = plt.figure()
ax = fig.add_subplot(111)
cursor = Cursor(ax, useblit=True, color='white', linewidth=0.5)
plt.subplots_adjust(bottom=0.2)
done_option = plt.axes([0.1, 0.05, 0.1, 0.075])
ax.set_xlim((-2, 2))
ax.set_ylim((-2, 2))
ax.set_title('click to build line segments')
ax.set_xlim((-1,2))
ax.set_ylim((0,5))

#cnv = Canvas(fig,ax,save_option,add_option,edit_option,plot_option,run_option,edit_shape_option,delete_option,split_option,X_length,Y_length,X_resolution)
linebuilder = Split_axis(fig,ax,done_option)
done = Button(done_option, 'done',color="black")
done.on_clicked(linebuilder.done)
plt.show()
if plt.get_fignums():
    print "asgdh"
else:
    print "figure is closed"
    print linebuilder.tx
"""