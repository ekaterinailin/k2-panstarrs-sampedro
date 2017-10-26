from astropy.coordinates import ICRS
from astropy import units as u
from astropy.coordinates.angles import Angle
from astropy.coordinates import match_coordinates_sky
import numpy as np
import os
import sys

class OpenCluster():

	'''
	Open cluster class.
	'''

	def __init__(self, name, sampedro_name, radius, age):
		
		self.name=name
		self.sampedro_name=sampedro_name
		self.radius=radius
		self.age=age
		self.sampedro_n0=[]
		self.sampedro_n1=[]
		self.sampedro_n2=[]
		self.sampedro_n3=[]
		self.PS=[]
		self.K2=[]
		self.K2MASS=[]
	
	def loadcatalog(self,path,ext,usecols,delimiter, skip_header, dtype):

		catfile=self.sampedro_name + ext

		if find(catfile,path):
		
			cat=np.genfromtxt(catfile ,delimiter=delimiter,skip_header=skip_header,usecols=usecols,dtype=dtype)
			ra=cat['ra'].tolist()
			dec=cat['dec'].tolist()
			ID=cat['ID'].tolist()
			
			if ('M1' and 'M2' and 'M3') in cat.dtype.names: 
				
				m1=cat['M1'].tolist()
				m2=cat['M2'].tolist()
				m3=cat['M3'].tolist()
				
				return [ID,ra,dec,m1,m2,m3]
			
		
			elif ('2MASS') in cat.dtype.names: 

				Twomass=cat['2MASS'].tolist()
			
				return [ID,ra,dec,Twomass]
			
			else:
				
				return [ID,ra,dec]

		else: 
			
			print('A ' + self.name + ' file is missing! Abort! Abort!')
			exit()
			
			return
		
	def loadcatalogs(self):
		
		path=sys.argv[1]
		os.chdir(path)
		
		print('Currently working in \"' + os.getcwd()+ '\"')

		self.sampedro_n0=self.loadcatalog(path,'_Sampedro_cluster_members_query.csv',(2,3,5,57,58,59),',',1 ,[('ID','U10'),('ra','f8'),('dec','f8'),('M1','b'),('M2','b'),('M3','b')])
		#workaround for weird RA in Sampedro...:
		self.sampedro_n0[1]=[x - 270.0 for x in self.sampedro_n0[1]] #weird coordinates...
		self.PS=self.loadcatalog(path, '_panstarrs_search.txt',(1,2,3),'\t',2 ,[('ID','i8'),('ra','U12'),('dec','U12')])
		self.K2=self.loadcatalog(path, '_k2_search.txt',(0,4,5,21),'\t',2 ,[('ID','i8'),('ra','U12'),('dec','U12'),('2MASS','U12')])
		
		#create the subset of K2 LCs where also 2MASS IDs exist
		#first of all copy the entire K2 list
		self.K2MASS=list(self.K2) 
		#empty string incidates there's no 2MASS ID
		no2mass = '' 
		#find all the rows where this string is found in the 2MASS ID column
		no2mass_idx=[i for i, x in enumerate(self.K2MASS[3]) if x==no2mass]
		
		#remove all rows with objects with no 2MASS ID
		for _ in range(len(K2.MASS)):
			self.K2MASS[_]=[i for j, i in enumerate(self.K2MASS[_]) if j not in no2mass_idx]

		return

	def refinesampedro(self):
	
		if self.sampedro_n0!=[]:
			l=self.sampedro_n0
			self.sampedro_n1=list(l)
			dell=[]
			for _ in range(len(l[0])):

				if (l[3][_] + l[4][_] + l[5][_]) < 1:
					dell.append(_)

			for _ in range(len(l)):
	
				self.sampedro_n1[_] = [i for j, i in enumerate(self.sampedro_n1[_]) if j not in dell]
		
		if self.sampedro_n1!=[]:
			l=self.sampedro_n1
			self.sampedro_n2=list(l)
			
			dell=[]
			for _ in range(len(l[0])):

				if (l[3][_] + l[4][_] + l[5][_]) < 2:
					dell.append(_)

			for _ in range(len(l)):
	
				self.sampedro_n2[_] = [i for j, i in enumerate(self.sampedro_n2[_]) if j not in dell]
			
		if self.sampedro_n2!=[]:
			self.sampedro_n3=list(self.sampedro_n2)
			dell=[]
			for _ in range(len(self.sampedro_n2[0])):

				if (self.sampedro_n2[3][_] + self.sampedro_n2[4][_] + self.sampedro_n2[5][_]) < 3:
					dell.append(_)

		
			if len(dell)==len(self.sampedro_n3[0]):
				print('Three methods do never seem to match?')
				self.sampedro_n3=[[0],[0],[0],[0],[0],[0]]
			
			else:
				for _ in range(len(self.sampedro_n2)):
	
					self.sampedro_n3[_] = [i for j, i in enumerate(self.sampedro_n3[_]) if j not in dell]
		
		
		print(len(self.sampedro_n0[5]))
		print(len(self.sampedro_n1[5]))
		print(len(self.sampedro_n2[5]))
		print(len(self.sampedro_n3[5]))

		return

	def sampedro_match(self, n, dist='0h0m3s', cat='Pan-STARRS'):

		'''
		
		Matches Sampedro catalog with a number (n) of membership approvals with catalog (cat) entries in a radius (dist) given in 'XhYmZs'.

		Output:

		match_idx: list - indices of objects in Pan-STARRS that match with object positions in Sampedro cluster list 

		'''
		
		sampedro=[self.sampedro_n0,self.sampedro_n1,self.sampedro_n2,self.sampedro_n3]
		
		if cat=='Pan-STARRS':
			l=self.PS	
		elif cat == 'K2MASS':
			l=self.K2MASS
		else:
			print('This is not a valid catalog chiffre. I\'ll use Pan-STARRS.')
			l=self.PS
		
		ra2=Angle(l[1],unit=u.deg)
		dec2=Angle(l[2],unit=u.deg)
			
		ra3=Angle(sampedro[n][1],u.deg) 
		dec3=Angle(sampedro[n][2],u.deg)

		#assume ra1/dec1 and ra/dec2 are arrays loaded from some file
		
		c = ICRS(ra=ra3, dec=dec3)#c=cluster stars from Sampedro
		catalog = ICRS(ra=ra2, dec=dec2)#catalog=panstarrs has all the parameters

		#match Sampedro cluster members with Pan-STARRS
		
		idx, d2d, d3d = match_coordinates_sky(c, catalog) #idx are indices into catalog that are the closest objects to each of the coordinates in c
		match_idx=are_within_bounds(d2d,'0h0m0s', dist)
		
		print(len(match_idx))
	
		return match_idx


#General purpose funcs:

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def are_within_bounds(arr, min_angle, max_angle):
	l=[]
	for _ in range(len(arr)):	
		if arr[_].is_within_bounds(min_angle, max_angle): l.append(_)
	return l


#Wrapping it:

x=OpenCluster('Ruprecht 147','Ruprecht_147', 30, 2.5)
x.loadcatalogs()
x.refinesampedro()
# print(x.sampedro_n0[1][:10])
# x.PS_sampedro_match(0)
# x.PS_sampedro_match(0,dist='0h0m5s')
# print(x.sampedro_n1[1][:10])
# x.PS_sampedro_match(1)
# x.PS_sampedro_match(1,dist='0h0m5s')
# print(x.sampedro_n2[1][:10])
# x.PS_sampedro_match(2)
# x.PS_sampedro_match(2,dist='0h0m5s')
# print(x.sampedro_n3[1][:10])
# x.PS_sampedro_match(3)
# x.PS_sampedro_match(3,dist='0h0m5s')

print(x.sampedro_n0[1][:10])
x.PS_sampedro_match(0)
x.PS_sampedro_match(0,dist='0h0m5s', cat='K2MASS')
print(x.sampedro_n1[1][:10])
x.PS_sampedro_match(1)
x.PS_sampedro_match(1,dist='0h0m5s', cat='K2MASS')
print(x.sampedro_n2[1][:10])
x.PS_sampedro_match(2)
x.PS_sampedro_match(2,dist='0h0m5s', cat='K2MASS')
print(x.sampedro_n3[1][:10])
x.PS_sampedro_match(3)
x.PS_sampedro_match(3,dist='0h0m5s', cat='K2MASS')


# print(type(K2['2MASS']))
# #refine K2 list by 2MASS condition
# TwoMASS=K2['2MASS']
# print(TwoMASS[:10])
# x=np.where(TwoMASS=='')
# print(len(TwoMASS))
# print(len(x[0]))
# x=np.where(TwoMASS!='')
# print(len(x[0]))
# 
# 
# #convert ra/dec into ICRS format
# ra1=Angle(K2['ra'],unit=u.deg)
# dec1=Angle(K2['dec'],unit=u.deg)
# ra2=Angle(PS['ra'],unit=u.deg)
# dec2=Angle(PS['dec'],unit=u.deg)
# ra3=Angle(SOC['ra']-270.0,u.deg) #weird coordinates...
# dec3=Angle(SOC['dec'],u.deg)
# 
# #assume ra1/dec1 and ra/dec2 are arrays loaded from some file
# c = ICRS(ra=ra3, dec=dec3)#c=cluster stars from Sampedro
# catalog = ICRS(ra=ra2, dec=dec2)#catalog=panstarrs has all the parameters
# k2_list=ICRS(ra=ra1, dec=dec1)
# 
# 
# #match Sampedro cluster members with Pan-STARRS
# idx, d2d, d3d = match_coordinates_sky(c, catalog) #idx are indices into catalog that are the closest objects to each of the coordinates in c
# #match Sampedro cluster members with K2 conditional on 2MASS
# idx2, d2d2, d3d2 = match_coordinates_sky(c, k2_list) #idx are indices into catalog that are the closest objects to each of the coordinates in c


#Find the closest object within 3 arcsec (following Huber et al (2016) for their matching procedure
#print(are_within_bounds(d2d,'0h0m0s', '0h0m30s'))
#print(are_within_bounds(d2d,'0h0m0s', '0h0m20s'))
#print(are_within_bounds(d2d,'0h0m0s', '0h0m10s'))
#print(len(are_within_bounds(d2d,'0h0m0s', '0h0m05s')))
#print(len(are_within_bounds(d2d,'0h0m0s', '0h0m03s')))
#print(are_within_bounds(d2d2,'0h0m0s', '0h0m30s'))
#print(are_within_bounds(d2d2,'0h0m0s', '0h0m20s'))
#print(are_within_bounds(d2d2,'0h0m0s', '0h0m10s'))
#print(len(are_within_bounds(d2d2,'0h0m0s', '0h0m05s')))
#print(len(are_within_bounds(d2d2,'0h0m0s', '0h0m03s')))






