from astropy.coordinates import ICRS
from astropy import units as u
from astropy.coordinates.angles import Angle
from astropy.coordinates import match_coordinates_sky
import random as rnd
import numpy as np
import os
import sys

def mock_wrap(size=1,dist=3):

	'''

	Main wrapper for catalog mockup creation.

	'''
	#The size of the mockup cat should be div'able into fractions down to 5%:
	size=int(size)*20
	#Mock a K2 catalog list of 1*size including 0.25*size 2MASS IDs
	#and hand over 0.1*size matches from K2MASS to Sampedro:
	k2,twomass_sampedro_match_idx=mock_k2(size)
	IDprint(k2, catname='K2')	
	#Mock a Pan-STARRS list of 2*size with 0.5*size matching K2
	ps=mock_panstarrs(size,k2,dist)
	IDprint(ps, catname='Pan-STARRS')
	#Mock a Sampedro list of 2*size with 0.1*size matching K2MASS, 	
	smp=mock_sampedro(size,k2,ps,twomass_sampedro_match_idx,dist)
	IDprint(smp, catname='Sampedro')
	for item in k2:
		print(item['ra'])

def randangle(val=0,stype='ra'):
	if val==0:
		r=rnd.gauss(15.,2.)
	else: 
		r=val	
	ra=Angle(r,unit=u.deg)
	if stype=='ra':
		ra=ra.to_string(sep='  ',pad=True)
	elif stype=='dec':
		ra=ra.to_string(sep='  ',pad=True,)
		ra=np.core.defchararray.add('+',ra)
	print(ra)
	return r, ra 


def mock_k2(size):
	
	'''

	Mocks a K2 catalog list including 25% fake 2MASS IDs.

	'''
	
	#tabseparated
	#28 columns (0,4,5,20) for ID, ra, dec, 2mass
	#example:
	#219622634	19.277982	-16.278071	19164073-1616411
	#first two rows occupied

	l=[]
	twomass_idx=[]
	for _ in range(size):
		if (_+1)%4==0: 
			tm='1'
			ID='K2MASS'
			twomass_idx.append(_)	
			r,ra=randangle()
			d,dec=randangle(stype='dec')		
		else: 
			tm=''
			ID='K2'
			r,ra=randangle()
			d,dec=randangle(stype='dec')
		l.append({"ID": ID, "ra": ra, "dec":dec, "twomass":tm,'instr_r':r,'instr_d':d})

	
	writecat(l,(0,4,5,20),name='mockup_k2_search.txt',tag='k2')
	#pick matches for Sampedro from those objects with 2MASS ID (40%)
	twomass_sampedro_match_idx=twomass_idx[:int(len(twomass_idx)*0.4)]

	return l,twomass_sampedro_match_idx

def mock_sampedro(size,k2,ps,twomass_sampedro_match_idx, dist):
	
	'''

	Mocks a Sampedro catalog list.

	'''
	
	#ID, ra, dec, M1, M2, M3
	l=[]
	
	#Fill in the matches from K2MASS

	for item in twomass_sampedro_match_idx:
		approx_ra=k2[item]['instr_r']+rnd.uniform(0.,dist/3600.)
		approx_dec=k2[item]['instr_d']+rnd.uniform(0.,dist/3600.)
		r,approx_ra=randangle(val=approx_ra)
		d,approx_dec=randangle(val=approx_dec,stype='dec')
		#l.append({"ID": 'K2MASS+Sampedro', "ra": approx_ra, "dec":approx_dec, "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})
		l.append({"ID": 'K2MASS+Sampedro', "ra": r, "dec":d, "M1":1,"M2":1,"M3":1,'instr_r':r,'instr_d':d})

	k2_match_start=len(l)

	#Remove matches from k2 list to avoid duplicates:

	k2 = [i for j, i in enumerate(k2) if j not in twomass_sampedro_match_idx]

	#Fill in matches from K2 and Panstarrs without necessary K2MASS:

	fillcat(l,ps,size,dist,frac=0.5,ID='K2/Pan-STARRS+Sampedro')

	#Fill the rest of the cat randomly

	k2_match_finish=len(l)	

	while len(l)<2*size:
		r,ra=randangle()
		d,dec=randangle(stype='dec')
		#l.append({"ID": 'Sampedro', "ra": rnd.gauss(25.,5.), "dec":rnd.gauss(25.,5.), "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})	
		l.append({"ID": 'Sampedro', "ra": r, "dec":d, "M1":1,"M2":1,"M3":1,'instr_r':r,'instr_d':d})
	
	writecat(l,(2,3,5,57,58,59),skiprows=1,name='mockup_Sampedro_cluster_members_query.csv',tag='smp')
	return l

def mock_panstarrs(size, k2,dist):
	
	'''

	Mocks a Pan-STARRS catalog list.

	'''

	#This is the catalog as a initialized list:

	l=[]

	#Fill 0.25*size with stars from K2

	fillcat(l,k2,size,dist,frac=0.25,ID='PS+K2')

	#Fill in the rest

	while len(l) < 2*size:
		r,ra=randangle()
		d,dec=randangle(stype='dec')
		l.append({"ID": 'PS', "ra": ra, "dec":dec,'instr_r':r,'instr_d':d})

	#Write the catalog into a file to use as debugging agent

	writecat(l,(1,2,3),name='mockup_panstarrs_search.txt',tag='ps')
	print('length of cat#2')
	print(l[1]['dec'])
	
	return l

def IDprint(cat, catname='Unknown'):

	'''

	Prints the IDs present in a mockup catalog.

	'''

	print(catname+':')
	for i in range(len(cat)): print(cat[i]['ID'])

	return 


def writecat(cat,usecols,skiprows=2,name='mockup_unknown.txt', tag='unknown'):
	
	'''

	Writes out calatogs into dummy files.

	'''
	#Make a copy of catalog for printing, so you can shuffle it a bit without interfering with the rest of the code:

	printcat=list(cat)	
	rnd.shuffle(printcat)

	#Create output file

	out=open('cats/'+name,'w')

	#Create a dummy entry to fill up remaining rows and columns

	dummy='dummy\t'	
	numcols=usecols[-1]+1

	#write the skipped rows

	dummymap=dummy*numcols
	for i in range(1,skiprows): dummymap+='\n'+dummy*numcols
	out.write(dummymap)

	#check specific catalog

	if tag=='k2': dicttag=['ID','ra','dec','twomass']
	elif tag=='ps': dicttag=['ID','ra','dec']
	elif tag=='smp': dicttag=['ID','ra','dec','M1','M2','M3']
	else: print('Use a valid tag.')

	#create line to write in file and write it

	for item in printcat:
		catline='\n'
		inkr=0
		for i in range(numcols):
			
			if i in usecols:
				catline+=str(item[dicttag[inkr]])+'\t'
				inkr+=1
			else: 
				catline+=dummy
		out.write(catline)
	out.close()
	return

	
def fillcat(l,cat,size,dist,frac=0.25,ID='ID'):

	'''

	Fills a catalog list _l_ with entries from another list _cat_ up to a certain fraction _frac_ of the initial _size_ measure and assigns these new entries _ID_. An uncertainty is added by shifting the coordinates by up to _dist_ arcmin.

	'''
	print('Length of cat: ')
	print(len(cat[1]))
	for i in range(int(frac*size)):
		print(cat[i]['instr_r'],cat[i]['instr_d'])		
		approx_ra=cat[i]['instr_r']+rnd.uniform(0.,dist/3600.)
		approx_dec=cat[i]['instr_d']+rnd.uniform(0.,dist/3600.)
		print(approx_ra,approx_dec)
		r,approx_ra=randangle(val=approx_ra)
		d,approx_dec=randangle(val=approx_dec,stype='dec')
		if ID=='K2/Pan-STARRS+Sampedro':
			#l.append({"ID": ID, "ra": approx_ra, "dec":approx_dec, "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})	
			l.append({"ID": ID, "ra": r, "dec":d, "M1":1,"M2":1,"M3":1,'instr_r':r,'instr_d':d})		
		elif ID=='PS+K2':
			l.append({"ID": ID, "ra": approx_ra, "dec":approx_dec,'instr_r':r,'instr_d':d})
		elif ID=='Sampedro':
			l.append({"ID": ID, "ra": r, "dec":d,'instr_r':r,'instr_d':d})

	
#WRAP!
mock_wrap()


