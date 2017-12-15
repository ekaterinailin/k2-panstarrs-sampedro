#https://matplotlib.org/2.1.0/users/event_handling.html

import numpy as np
import matplotlib.pyplot as plt
import sys
import datetime 

#--------------------------------------------------------------------------------------------------
#FUNCTIONS:
#--------------------------------------------------------------------------------------------------

def reset_pipeline():
	global count, myflare_start, myflare_stop, myflare_start_flux, myflare_stop_flux,comment,interrupted
	count=0#length of list instead?
	comment=''
	myflare_start,myflare_stop,myflare_start_flux,myflare_stop_flux=[],[],[],[]



def on_pick(event):  #event => matplotlib.backend_bases.PickEvent

	#Is there a more elegant solution instead of defining these variables as global?
	
	global count, myflare_start, myflare_stop, myflare_start_flux, myflare_stop_flux,comment,interrupted
	
	#What happens if I click on a data point:

	#Line is picked:

	thisline = event.artist #matplotlib.lines.Line2D => A line in matplotlib
	
	#The chosen data point is extracted:

	#First info for ALL data point on this line are extracted:

	xdata = thisline.get_xdata()
	ydata = thisline.get_ydata()

	#Then the specific event is chosen by index and the data are read out:

	ind = event.ind
	time = xdata[ind]
	flux = ydata[ind]
	
	#Sometimes the choice by clicking is ambiguous:

	if len(time)>1:
		print('ATTENTION!\nYou picked more than one data point at once. Try again.\nATTENTION\n')
	
	#If not, I can safely interact with the figure as I wish:

	else:
		print('onpick time:', time)

		# I only want to save flares that have START and END, so I pick events in pairs and count on the go:
		count+=1
		
		#If I only have one event, I add it to the starts of flares list:
		
		if count%2==1:

			myflare_start.append(time)
			myflare_start_flux.append(flux)

			#I also want to see the point I registered as start of a flare:

			ax.plot(myflare_start,myflare_start_flux,'o',alpha=0.8, lw=0.4,color='green')

		#If I collected two events, I can choose to save this pair as a flare and write the data into a file:

		elif count%2==0:

			myflare_stop.append(time)
			myflare_stop_flux.append(flux)

			#I also want to see the point I registered as start of a flare:

			ax.plot(myflare_stop,myflare_stop_flux,'o',alpha=0.8, lw=0.4,color='red')
			
			#Here I may trigger a save-and-proceed or remove-and-proceed event:

			print('Press \"y\" to confirm flare events. Press \"x\" to remove.')
	
	#This command shows the changes to the figure:

	plt.draw()

	return

#If I trigger a save-and-proceed or remove-and-proceed event I end up here:

def on_key(event):
	
	#Here I write in my pairs of flare starts and ends:

	myflares=open('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'_my_flares.txt','a')
	if interrupted==True:
		myflares.write(',,,,'+str(datetime.datetime.now())+'\n')
		interrupted=False
	#Again: is there a more elegant solution?

	global count, myflare_start, myflare_stop, myflare_start_flux, myflare_stop_flux,comment,interrupted, time	
	#So that you know:

	print('You pressed', event.key)
		
	#If I choose to save the marked pairs:

	if event.key=='y':
		
		print('The following events are added to the list:\n')
		for i in range(count//2):
			line=str(myflare_start[i][0])+','+ str(myflare_stop[i][0])
			print(line)
			print()
		try:
			comment=input('Click into command prompt to write and enter a comment on the selected events,\nthen click into the figure\nand/or press \"a\": \n')
		except RuntimeError:
			print('You failed.')
		print()
		
	#Elif I choose to discard and try again:

	elif event.key=='x':

		print('The following events are removed from the list:\n')
		for i in range(count//2):

			ax.plot(myflare_stop,myflare_stop_flux,'o',alpha=0.8, lw=0.4,color='blue')
			ax.plot(myflare_start,myflare_start_flux,'o',alpha=0.8, lw=0.4,color='blue')
			plt.draw()
			line=str(myflare_start[i][0])+','+ str(myflare_stop[i][0])+'\n'
			print(line)
		#In any case I set back the counter and empty the lists:
		reset_pipeline()
		myflares.close()
		
	#If I choose to comment or not I write out the following line:

	elif event.key=='a':
		for i in range(count//2):
			print('This is your comment: '+comment+'\n')
			line=str(myflare_start[i][0])+','+ str(myflare_stop[i][0])+','+str(time.index(myflare_start[i][0]))+','+str(time.index(myflare_stop[i][0]))+','+str(comment)
			print(line)
			myflares.write(line+'\n')
		
		#In any case I set back the counter and empty the lists:
		del comment
		reset_pipeline()
		myflares.close()

	return

def loaddata(cluster, objectid):
	
	#Read in data and convert to lists:
	time,flux_gap,error,flux_model=np.loadtxt('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'.txt',delimiter=',',unpack=True)
	istart, istop=np.loadtxt('share/cross_match/union/'+cluster+'_post_appa/'+objectid+'_flares.txt',delimiter=',',unpack=True,dtype=np.dtype(np.int16))

	return list(time),list(flux_gap),list(error),list(flux_model),list(istart), list(istop)


def plotLC(time,flux_gap,error,flux_model,istart, istop,objectid):

	global interrupted
	#Set up a plot with matplotlib:
	A=plt.figure()
	ax =plt.subplot(211)
	ax.plot(time, flux_gap,'o',alpha=0.8, lw=0.4, picker=5) #observation
	ax.errorbar(time, flux_gap,yerr=error,alpha=0.8, lw=1) #statistical error
	for g,start in enumerate(istart):
		xregion=time[istart[g]:istop[g]+1]
		yregion=flux_gap[istart[g]:istop[g]+1]
		ax.plot(xregion,yregion,color='red', lw=2) #detected flares in observations
		#Shade the region where a flare is detected by Appaloosa:
		ax.axvspan(xregion[0],xregion[-1],edgecolor='black', alpha=0.2,linewidth=2)

	ax.plot(time, flux_model, 'blue', lw=0.5) #model light curve from which the flare signatures deviate

	#Cosmetics

	plt.title('EPIC '+ objectid)
	plt.ylabel(r'Flux ($e^-$ sec$^{-1}$)')

	#Plot range specification

	xdur0 = min(time)
	xdur1 = max(time)
	xdurok = np.where((time >= xdur0) & (time <= xdur1))
	xdurok=xdurok[0] #xdurok is a tuple with a list: conversion needed
	plt.xlim(xdur0, xdur1) 
	plt.ylim(min([flux_gap[x] for x in xdurok]), max([flux_gap[x] for x in xdurok]))


	#Set up a residuals plot with matplotlib:

	bx=plt.subplot(212,sharex=ax)
	residuals=[x-y for x,y in zip(flux_gap,flux_model)]
	bx.plot(time, residuals,alpha=0.8, lw=1) #residuals
	plt.xlim(xdur0, xdur1)
	plt.ylim(max(-500,min([residuals[x] for x in xdurok])), max([500]+[residuals[x] for x in xdurok]))
	for g,start in enumerate(istart):
		xregion=time[istart[g]:istop[g]+1]
		yregion=flux_gap[istart[g]:istop[g]+1]
		#Shade the region where a flare is detected by Appaloosa:
		bx.axvspan(xregion[0],xregion[-1],edgecolor='black', alpha=0.2,linewidth=2)
	plt.xlabel('Time (BJD - 2454833 days)')
	plt.ylabel(r'Residual Flux ($e^-$ sec$^{-1}$)')

	return A,ax,bx

def wrap(cluster, objectid,interrupted=False):

	#Load and plot data:

	time,flux_gap,error,flux_model,istart, istop=loaddata(cluster, objectid)
	reset_pipeline()
	plot,ax,bx=plotLC(time,flux_gap,error,flux_model,istart, istop,objectid)
	global ax,time
	#Here is where I link the interactive figure to the code:
	#The FigureCanvas method mpl_connect() returns a connection id which is simply an integer.

	#Connect event with string 'key_press_event' to function on_key:

	cid = plot.canvas.mpl_connect('key_press_event', on_key) 

	#Connect event with string 'pick_event' to function on_pick:

	cid2=plot.canvas.mpl_connect('pick_event', on_pick)

	plt.show()

	return

#--------------------------------------------------------------------------------------------------
#EXECUTION:
#--------------------------------------------------------------------------------------------------

#Kepler EPIC ID for the light curve I want to look at:
cluster='Ruprecht_147'
#objectid='219396051'

#Do I want to analyse this LC or have I already done so?
import glob, os
working_directory = os.getcwd()
EPIC=open('share/cross_match/union/'+cluster+'_union_IDs_wo_prefix.txt')
os.chdir("share/cross_match/union/"+cluster+"_post_appa")
interrupted=False
for ID in EPIC:
	objectid=ID[:9]
	global objectid
	print(glob.glob(ID+"_my_flares.txt"))
	if glob.glob(objectid+"_my_flares.txt")==[]:
		os.chdir(working_directory)
		wrap(cluster, objectid)
	else:
		user=input('You have at least started analysing LC ' + objectid+'. Want to work on it anyway? (y/n)')
		if user=='y':
			os.chdir(working_directory)
			interrupted=True
			wrap(cluster, objectid, interrupted=interrupted)
		

