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

	def __init__(self, name, sampedro_name, radius=0., age=0.):
		
		#prompt and file specific cluster name
		self.name=name
		self.sampedro_name=sampedro_name
		#radius in acrmin
		self.radius=radius
		#age in Gyr
		self.age=age
		#subcats with different numbers of membership allocations
		self.sampedro_n0=[]
		self.sampedro_n1=[]
		self.sampedro_n2=[]
		self.sampedro_n3=[]
		#Pan-STARRS catalog
		self.PS=[]
		#K2 search for the cluster 
		self.K2=[]
		#
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
				#print([ID,ra,dec,m1,m2,m3])
				return [ID,ra,dec,m1,m2,m3]
			
		
			elif ('2MASS') in cat.dtype.names: 

				Twomass=cat['2MASS'].tolist()
				#print([ID,ra,dec,Twomass])
				return [ID,ra,dec,Twomass]
			
			else:
				#print([ID,ra,dec])
				return [ID,ra,dec]

		else: 
			
			print('A ' + self.name + ' file is missing! Abort! Abort!')
			exit()
			
			return
		
	def loadcatalogs(self):
		
		path=sys.argv[1]
		os.chdir(path)
		
		print('Currently working in \"' + os.getcwd()+ '\"')

		self.sampedro_n0=self.loadcatalog(path,'_Sampedro_cluster_members_query.csv',(2,3,5,57,58,59),'\t',1 ,[('ID','U10'),('ra','f8'),('dec','f8'),('M1','b'),('M2','b'),('M3','b')])
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
		for _ in range(len(self.K2MASS)):
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
			l=list(self.PS)	
		elif cat == 'K2MASS':
			l=list(self.K2MASS)
		else:
			print('This is not a valid catalog chiffre. I\'ll use Pan-STARRS.')
			l=self.PS
		
		ra2=Angle(l[1],unit=u.deg)
		dec2=Angle(l[2],unit=u.deg)
			
		ra3=Angle(sampedro[n][1],u.deg) 
		dec3=Angle(sampedro[n][2],u.deg)

		#assume ra1/dec1 and ra/dec2 are arrays loaded from some file
		
		c = ICRS(ra=ra3, dec=dec3)#cluster stars from Sampedro
		catalog = ICRS(ra=ra2, dec=dec2)#panstarrs or K2MASS has all the parameters

		#match Sampedro cluster members with Pan-STARRS
	#	print(len(c))
	#	print(len(catalog))
		idx, d2d, d3d = match_coordinates_sky(c,catalog) #idx are indices into catalog that are the closest objects to each of the coordinates in c
	#	print(idx)
	#	print(d2d)
		match_idx=are_within_bounds(idx,d2d,'0h0m0s', dist)
	#	print(match_idx)
		print('Sampedro_n' + str(n) + ' X ' + cat + ' X ' + dist + ': ' + str(len(match_idx)) + ' matching objects.')
	#	print(len(match_idx))
	
		return match_idx

	def second_order_match(self, n, distPS='0h0m3s', dist='0h0m3s'):

		'''
		
		Matching K2 and Pan-STARRS (conditional on Sampedro cluster membership with (n) shared assessments).
		
		'''

		lK2=list(self.K2)
		lPS=list(self.PS)

		match_idx_PS=self.sampedro_match(n, dist=distPS, cat='Pan-STARRS')

		idx_K2=list(range(len(lK2[0])))#is that one even necessary???
		#print('K2 object list has '+ str(len(idx_K2))+ ' LCs.')

		for _ in range(len(self.K2)):
			lK2[_]= [i for j, i in enumerate(lK2[_]) if j in idx_K2]

		for _ in range(len(self.PS)):
			lPS[_]= [i for j, i in enumerate(lPS[_]) if j in match_idx_PS]

	#	print('From the ' + self.name + ' members as assessed by Sampedro ' + str(len(lPS[0])) + ' have an Pan-STARRS ID.')
		ra2=Angle(lK2[1],unit=u.deg)
		dec2=Angle(lK2[2],unit=u.deg)
			
		ra3=Angle(lPS[1],u.deg) 
		dec3=Angle(lPS[2],u.deg)

		#assume ra1/dec1 and ra/dec2 are arrays loaded from some file
		
		c = ICRS(ra=ra3, dec=dec3)#c=cluster stars from Pan-STARRS
		catalog = ICRS(ra=ra2, dec=dec2)#catalog=K2 is where I need LCs to exist

		#match Sampedro cluster members with Pan-STARRS
		
		idx, d2d, d3d = match_coordinates_sky(c, catalog,nthneighbor=1) #idx are indices into catalog that are the closest objects to each of the coordinates in c
		second_match_idx=are_within_bounds(idx, d2d,'0h0m0s', dist)
		print('Sampedro_n' + str(n) + ' X  Pan-STARRS X K2 X ' + dist + ': ' + str(len(second_match_idx)) + ' matching objects.')
		
		return second_match_idx



#General purpose funcs:

def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def are_within_bounds(idx,d2d, min_angle, max_angle):
	l=[]
	d=[]
	for _ in range(len(d2d)):	
		if d2d[_].is_within_bounds(min_angle, max_angle): 
			l.append(idx[_])
			d.append(d2d[_])
	
	n=[]
	lpop=list(l)
	dpop=list(d)
	while len(lpop)>0:
		#find same indices
		m=[i for i, j in enumerate(lpop) if j == lpop[0]]
	#	print('m:')
	#	print(len(m))
		if len(m)>1:
	#		print(m)
			for k in range(len(m)-1):
				if dpop[m[k]].is_within_bounds(dpop[m[k+1]]): 
					minidx=lpop[m[k]]
				else: 
					minidx=lpop[m[k+1]]
		else: minidx=lpop[m[0]]
		n.append(minidx)
		lpop=[i for j, i in enumerate(lpop) if j not in m]
		dpop=[i for j, i in enumerate(dpop) if j not in m]
		# print('lpop:')
		# print(lpop)
		# print('n:')
		# print(n)
	

	return n


#Wrapping it:

x=OpenCluster('mockup','mockup', radius=30, age=2.5)
x.loadcatalogs()
x.refinesampedro()



for i in range(1,4):
	#10
	x.sampedro_match(i,dist='0h0m5s',cat='Pan-STARRS')
	#11
	x.sampedro_match(i,dist='0h0m3s',cat='Pan-STARRS')
for i in range(1,4):
	#16,18,20
	x.sampedro_match(i,dist='0h0m5s',cat='K2MASS')
	#17,19,21     
	x.sampedro_match(i,dist='0h0m3s',cat='K2MASS')
for i in range(1,4):
	#22,26,30
	x.second_order_match(i,distPS='0h0m5s',dist='0h0m5s')
	#23,27,31
	x.second_order_match(i,distPS='0h0m5s',dist='0h0m3s')
	#24,28,32
	x.second_order_match(i,distPS='0h0m3s',dist='0h0m5s')
	#25,29,33
	x.second_order_match(i,distPS='0h0m3s',dist='0h0m3s')

# for i in range(1,4):
#     #34,35,36
#     x.second_order_match(i,distPS='0h0m10s',dist='0h0m10s')
#     #34,35,36
#     x.second_order_match(i,distPS='0h0m5s',dist='0h0m5s')
# 

# print('Here\'s the best choice: ')
# A=x.sampedro_match(1,dist='0h0m3s',cat='K2MASS')
# results=[]
# # #print(A)
# results=[i for j, i in enumerate(x.K2MASS[0]) if j in A]
# print(results)
#print(len(results))
