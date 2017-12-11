import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('click on points')

line, = ax.plot(np.random.rand(100), 'o', picker=5)  # 5 points tolerance

def onpick(event):
	thisline = event.artist
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)

fig.canvas.mpl_connect('pick_event', onpick)
import matplotlib
print(matplotlib.get_backend())
plt.show()


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

