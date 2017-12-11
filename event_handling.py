import numpy as np
import matplotlib.pyplot as plt

objectid='205071984'
time,flux_gap,error,flux_model=np.loadtxt('share/cross_match/union/M44_post_appa/'+objectid+'.txt',delimiter=',',unpack=True)
istart, istop=np.loadtxt('share/cross_match/union/M44_post_appa/'+objectid+'_flares.txt',delimiter=',',unpack=True,dtype=np.dtype(np.int16))
time,flux_gap,error,flux_model,istart, istop=list(time),list(flux_gap),list(error),list(flux_model),list(istart), list(istop)

A = plt.figure(figsize=(5,1))
ax = A.add_subplot(111)
ax.plot(time, flux_gap,'o',alpha=0.8, lw=0.4, picker=5)#,ecolor='g')
#, yerr=error
for g,start in enumerate(istart):
	plt.plot(time[istart[g]:istop[g]+1],flux_gap[istart[g]:istop[g]+1],color='red', lw=1)

plt.plot(time, flux_model, 'blue', lw=0.5)
		
plt.title('EPIC '+ objectid, fontsize=8)
plt.xlabel('Time (BJD - 2454833 days)', fontsize=8)
plt.ylabel(r'Flux ($e^-$ sec$^{-1}$)', fontsize=8)

plt.xticks(fontsize=8, rotation=0)
plt.yticks(fontsize=8, rotation=0)
xdur = max(time) - min(time)
xdur0 = min(time)# + xdur/2.
xdur1 = max(time)# + xdur/1.5
xdurok = np.where((time >= xdur0) & (time <= xdur1))
xdurok=xdurok[0]

plt.xlim(xdur0, xdur1) # only plot a chunk of the data
print(xdurok)
plt.ylim(min([flux_gap[x] for x in xdurok]), max([flux_gap[x] for x in xdurok]))
#plt.savefig('share/cross_match/union/M44_post_appa/'+objectid + '_lightcurve.png', dpi=300, bbox_inches='tight', pad_inches=0.02)

def onpick(event):
	print(type(event))
	thisline = event.artist
	print(type(thisline))
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)

A.canvas.mpl_connect('pick_event', onpick)
plt.show()

