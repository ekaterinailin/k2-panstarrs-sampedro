import numpy as np
import matplotlib.pyplot as plt

#Kepler EPIC ID for the light curve I want to look at:
cluster='M44'
objectid='205071984'

#Read in data and convert to lists:

time,flux_gap,error,flux_model=np.loadtxt('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'.txt',delimiter=',',unpack=True)
istart, istop=np.loadtxt('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'_flares.txt',delimiter=',',unpack=True,dtype=np.dtype(np.int16))
time,flux_gap,error,flux_model,istart, istop=list(time),list(flux_gap),list(error),list(flux_model),list(istart), list(istop)


#Set up a plot with matplotlib:

#A = plt.figure(figsize=(5,1))
#ax = A.add_subplot(111)
A, ax =plt.subplots()
bx=ax.plot(time, flux_gap,'o',alpha=0.8, lw=0.4, picker=5) #observation
plt.errorbar(time, flux_gap,yerr=error,alpha=0.8, lw=1) #statistical error
for g,start in enumerate(istart):
	xregion=time[istart[g]:istop[g]+1]
	yregion=flux_gap[istart[g]:istop[g]+1]
	plt.plot(xregion,yregion,color='red', lw=2) #detected flares in observations
	plt.axvspan(xregion[0],xregion[-1],edgecolor='black', alpha=0.2,linewidth=2)

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


count=0
myflare_start,myflare_stop,myflare_start_flux,myflare_stop_flux=[],[],[],[]


def on_pick(event):  #event => matplotlib.backend_bases.PickEvent
	global count, myflare_start, myflare_stop, myflare_start_flux, myflare_stop_flux
	thisline = event.artist #matplotlib.lines.Line2D
	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()
	ind = event.ind
	time = xdata[ind]
	flux = ydata[ind]
	if len(time)>1:
		print('ATTENTION!\nYou picked more than one data point at once. Try again.\nATTENTION\n')
	else:
		print('onpick time:', time)
		count+=1

		if count%2==1:
			myflare_start.append(time)
			myflare_start_flux.append(flux)
			plt.plot(myflare_start,myflare_start_flux,'o',alpha=0.8, lw=0.4,color='green')


		elif count%2==0:
			myflare_stop.append(time)
			myflare_stop_flux.append(flux)
			plt.plot(myflare_stop,myflare_stop_flux,'o',alpha=0.8, lw=0.4,color='red')
			print('Press \"enter\" to confirm flare events. Press \"x\" to remove.')
		
	plt.draw()

	return

def on_key(event):
	global count, myflare_start, myflare_stop, myflare_start_flux, myflare_stop_flux
	myflares=open('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'_my_flares.txt','a')
	print('You pressed', event.key)
	
	if event.key=='enter':
	
		print('The following events are added to the list:\n')
		for i in range(count//2):
			line=str(myflare_start[i][0])+','+ str(myflare_stop[i][0])
			print(line)
			print()
			myflares.write(line+'\n')

	elif event.key=='x':

		print('The following events are removed from the list:\n')
		for i in range(count//2):

			plt.plot(myflare_stop,myflare_stop_flux,'o',alpha=0.8, lw=0.4,color='blue')
			plt.plot(myflare_start,myflare_start_flux,'o',alpha=0.8, lw=0.4,color='blue')
			line=str(myflare_start[i][0])+','+ str(myflare_stop[i][0])+'\n'
			print(line)
	
	count=0
	myflare_start,myflare_stop,myflare_start_flux,myflare_stop_flux=[],[],[],[]


			
	myflares.close()
	return


#Connect event with string 'key_press_event' to function on_key:
cid = A.canvas.mpl_connect('key_press_event', on_key) 
#Connect event with string 'pick_event' to function on_pick:
cid2=A.canvas.mpl_connect('pick_event', on_pick)#The FigureCanvas method mpl_connect() returns a connection id which is simply an integer.
plt.show()

