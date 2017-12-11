import numpy as np
import matplotlib.pyplot as plt

#Kepler EPIC ID for the light curve I want to look at:

objectid='205071984'

#Read in data and convert to lists:

time,flux_gap,error,flux_model=np.loadtxt('share/cross_match/union/M44_post_appa/'+objectid+'.txt',delimiter=',',unpack=True)
istart, istop=np.loadtxt('share/cross_match/union/M44_post_appa/'+objectid+'_flares.txt',delimiter=',',unpack=True,dtype=np.dtype(np.int16))
time,flux_gap,error,flux_model,istart, istop=list(time),list(flux_gap),list(error),list(flux_model),list(istart), list(istop)


#Set up a plot with matplotlib:

A = plt.figure(figsize=(5,1))
ax = A.add_subplot(111)
ax.plot(time, flux_gap,'o',alpha=0.8, lw=0.4, picker=5) #observation
plt.errorbar(time, flux_gap,yerr=error,alpha=0.8, lw=1) #statistical error
for g,start in enumerate(istart):
	plt.plot(time[istart[g]:istop[g]+1],flux_gap[istart[g]:istop[g]+1],color='red', lw=1) #detected flares in observations

plt.plot(time, flux_model, 'blue', lw=0.5) #model light curve from which the flare signatures deviate
		
#Cosmetics

plt.title('EPIC '+ objectid)
plt.xlabel('Time (BJD - 2454833 days)')
plt.ylabel(r'Flux ($e^-$ sec$^{-1}$)')

#Plot range specification

xdur0 = min(time)
xdur1 = max(time)
xdurok = np.where((time >= xdur0) & (time <= xdur1))
xdurok=xdurok[0] #xdurok is a tuple with a list: conversion needed
plt.xlim(xdur0, xdur1) 
plt.ylim(min([flux_gap[x] for x in xdurok]), max([flux_gap[x] for x in xdurok]))

def onpick(event):  #event => matplotlib.backend_bases.PickEvent
	thisline = event.artist #matplotlib.lines.Line2D
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	points = tuple(zip(xdata[ind], ydata[ind]))
	print('onpick points:', points)
	return points

f=A.canvas.mpl_connect('pick_event', onpick)
plt.show()

