import numpy as np
import matplotlib.pyplot as plt

objectid='205071984'
item,flux_gap,error,flux_model=np.loadtxt(objectid+'.txt',delimiter=',',unpack=True)
istart, istop=np.loadtxt(objectid+'_flares.txt',delimiter=',',unpack=True)
plt.subplots(figsize=(5,1))
plt.errorbar(time, flux_gap, yerr=error, color='k',alpha=0.8, lw=0.4,ecolor='g')

for g in range(len(istart)):
    plt.plot(time[istart[g]:istop[g]+1],flux_gap[istart[g]:istop[g]+1],color='red', lw=0.4)

plt.plot(time, flux_model, 'blue', lw=0.5)
        
plt.title('EPIC '+ objectid, fontsize=8)
plt.xlabel('Time (BJD - 2454833 days)', fontsize=8)
plt.ylabel(r'Flux ($e^-$ sec$^{-1}$)', fontsize=8)

plt.xticks(fontsize=8, rotation=0)
plt.yticks(fontsize=8, rotation=0)
xdur = np.nanmax(time) - np.nanmin(time)
#xdur0 = np.nanmin(time)# + xdur/2.
#xdur1 = np.nanmax(time)# + xdur/1.5
#xdurok = np.where((time >= xdur0) & (time <= xdur1))
plt.xlim(xdur0, xdur1) # only plot a chunk of the data
plt.ylim(np.nanmin(flux_gap[xdurok]), np.nanmax(flux_gap[xdurok]))
plt.savefig(file + '_lightcurve.png', dpi=300, bbox_inches='tight', pad_inches=0.02)

#def onpick(event):
	#thisline = event.artist
	#xdata = thisline.get_xdata()
	#ydata = thisline.get_ydata()
	#ind = event.ind
	#points = tuple(zip(xdata[ind], ydata[ind]))
	#print('onpick points:', points)

#fig.canvas.mpl_connect('pick_event', onpick)
#import matplotlib
#print(matplotlib.get_backend())
#plt.show()


# 
# import gtk
# 
# from matplotlib.figure import Figure
# from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
# from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
# 
# win = gtk.Window()
# win.connect("destroy", lambda x: gtk.main_quit())
# win.set_default_size(400,300)
# win.set_title("Embedding in GTK")
# 
# vbox = gtk.VBox()
# win.add(vbox)
# 
# fig = Figure(figsize=(5,4), dpi=100)
# ax = fig.add_subplot(111)
# ax.plot([1,2,3])
# 
# canvas = FigureCanvas(fig)  # a gtk.DrawingArea
# vbox.pack_start(canvas)
# toolbar = NavigationToolbar(canvas, win)
# vbox.pack_start(toolbar, False, False)
# 
# win.show_all()
# gtk.main()
# 

