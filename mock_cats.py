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
	size=int(size)*20
	#Mock a K2 catalog list of 1*size including 0.25*size 2MASS IDs
	#Hand over 0.1*size matches from K2MASS to Sampedro
	k2,twomass_sampedro_match_idx=mock_k2(size)
	print('K2:')
	for i in range(len(k2)): print(k2[i]['ID'])
	#Mock a Pan-STARRS list of 2*size with 0.5*size matching K2
	ps=mock_panstarrs(size,k2,dist)
	print('PS:')	
	for i in range(len(ps)): print(ps[i]['ID'])
	#Mock a Sampedro list of 2*size with 0.1*size matching K2MASS, 	
	smp, match=mock_sampedro(size,k2,ps,twomass_sampedro_match_idx,dist)
	print('SMP:')	
	for i in range(len(smp)): print(smp[i]['ID'])
	#write catalogs to .txt or .csv files


def mock_k2(size):
	
	'''

	Mocks a K2 catalog list including 25% fake 2MASS IDs.

	'''
	
	#tabseparated
	#28 columns 0,4,5,21 for ID, ra, dec, 2mass
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
		else: 
			tm=''
			ID='K2'
		l.append({"ID": ID, "ra": rnd.gauss(20.,5.), "dec":rnd.gauss(20.,5.), "twomass":tm})
		#print(l[_]['ID'])

	out=open('mockup_k2_search.txt','w')
	dummy='dummy\t'	
	dummymap=(dummy*28+'\n'+ dummy*28)
	out.write(dummymap)
	for item in l:
		out.write('\n'+str(item['ID'])+'\t'+dummy*3+str(item['ra'])+'\t'+str(item['dec'])+'\t' + dummy*15+str(item['twomass'])+'\t'+dummy*6)
	out.close()	
	
	#pick matches for Sampedro from those objects with 2MASS ID (40%)
	twomass_sampedro_match_idx=twomass_idx[:int(len(twomass_idx)*0.4)]
	#print(twomass_idx)
	#print(twomass_sampedro_match_idx)

	return l,twomass_sampedro_match_idx

def mock_sampedro(size,k2,ps,twomass_sampedro_match_idx, dist):
	
	'''

	Mocks a Sampedro catalog list.

	'''
	#ID, ra, dec, M1, M2, M3
	l=[]
	
	#Fill in the matches from K2MASS
	for item in twomass_sampedro_match_idx:
		approx_ra=float(k2[item]['ra'])+rnd.uniform(0.,dist/3600.)
		approx_dec=float(k2[item]['dec'])+rnd.uniform(0.,dist/3600.)
		l.append({"ID": 'K2MASS+Sampedro', "ra": approx_ra, "dec":approx_dec, "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})
		#print(k2[item])
	k2_match_start=len(l)
	#Remove matches from k2 list to avoid duplicates:
	k2 = [i for j, i in enumerate(k2) if j not in twomass_sampedro_match_idx]
	#Fill in matches from K2 and Panstarrs without necessary K2MASS:
	for i in range(int(0.5*size)):
		approx_ra=float(ps[i]['ra'])+rnd.uniform(0.,dist/3600.)
		approx_dec=float(ps[i]['dec'])+rnd.uniform(0.,dist/3600.)
		l.append({"ID": ps[i]['ID']+'+Sampedro', "ra": approx_ra, "dec":approx_dec, "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})
	#Fill the rest of the cat randomly
	k2_match_finish=len(l)	
	while len(l)<2*size:
		l.append({"ID": 'Sampedro', "ra": rnd.gauss(25.,5.), "dec":rnd.gauss(25.,5.), "M1":np.random.choice((0,1)),"M2":np.random.choice((0,1)),"M3":np.random.choice((0,1))})	

	sampedro_k2_match_idx=l[k2_match_start:k2_match_finish]
	print(len(sampedro_k2_match_idx))

	out=open('mockup_Sampedro_cluster_members_query.csv','w')
	#(2,3,5,57,58,59)
	dummy='dummy\t'	
	out.write(dummy*59)
	for item in l:
		out.write('\n'+dummy*2+str(item['ID'])+'\t'+str(item['ra'])+'\t'+dummy+str(item['dec'])+'\t' + dummy*51+str(item['M1'])+'\t'+str(item['M2'])+'\t'+str(item['M3']))
	out.close()	

	return l, sampedro_k2_match_idx

def mock_panstarrs(size, k2,dist):
	
	'''

	Mocks a Pan-STARRS catalog list.

	'''
	l=[]
	#Fill 0.25*size with stars from K2
	for i in range(int(0.25*size)):
		approx_ra=float(k2[i]['ra'])+rnd.uniform(0.,dist/3600.)
		approx_dec=float(k2[i]['dec'])+rnd.uniform(0.,dist/3600.)		
		l.append({"ID": 'PS+K2', "ra": approx_ra, "dec":approx_dec})
	#Fill in the rest
	while len(l) < 2*size:
		
		l.append({"ID": 'PS', "ra": rnd.gauss(20.,5.), "dec":rnd.gauss(20.,5.)})
		#print(l[_]['ID'])

	out=open('mockup_panstarrs_search.txt','w')
	dummy='dummy\t'	
	dummymap=(dummy*4+'\n'+dummy*4)
	out.write(dummymap)
	for item in l:
		out.write('\n'+dummy+str(item['ID'])+'\t'+str(item['ra'])+'\t'+str(item['dec']))
	out.close()	
	
	return l
#WRAP!
mock_wrap()

