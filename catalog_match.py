#Cleaning to-dos:
#clean up mocktest()
#write docstrings
#remove confusing comments


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
		#K2 EVEREST search for the cluster 
		self.K2=[]
		#subcat of self.K2 which has 2MASS IDs
		self.K2MASS=[]

	#----------------------------------------------------------------------------------------------
	
	def loadcatalog(self,path,ext,usecols,delimiter, skip_header, dtype,debug=False):

		catfile=self.sampedro_name + ext

		if find(catfile,path):
		
			cat=np.genfromtxt('cats/'+catfile ,delimiter=delimiter,skip_header=skip_header,usecols=usecols,dtype=dtype)
			ra=cat['ra'].tolist()
			dec=cat['dec'].tolist()
			ID=cat['ID'].tolist()
			
			if ('M1' and 'M2' and 'M3') in cat.dtype.names: 
				m1=cat['M1'].tolist()
				m2=cat['M2'].tolist()
				m3=cat['M3'].tolist()
				ra=Angle(ra,unit=u.deg)
				dec=Angle(dec, unit=u.deg)
				ra=ra.to_string(sep='dms',pad=True)
				dec=dec.to_string(sep='dms',pad=True,)
				return [ID,ra,dec,m1,m2,m3]
		
			elif ('2MASS') in cat.dtype.names:
				for item in ra:
					ra[ra.index(item)]='+'+item[:2]+'h'+item[3:5]+'m'+item[6:]+'s'
				for item in dec:
					dec[dec.index(item)]=item[:3]+'d'+item[4:6]+'m'+item[7:]+'s'
				Twomass=cat['2MASS'].tolist()
				return [ID,ra,dec,Twomass]
			
			else:

				#for item in ra:
				#	ra[ra.index(item)]='+'+item[:2]+'h'+item[3:5]+'m'+item[6:]+'s'
				#for item in dec:
				#	dec[dec.index(item)]=item[:3]+'d'+item[4:6]+'m'+item[7:]+'s'
				#print([ID,ra,dec])
				#ra=[float(x) for x in ra]
				#dec=[float(x) for x in dec]
				ra=Angle(ra,unit=u.deg)
				dec=Angle(dec, unit=u.deg)
				ra=ra.to_string(sep='dms',pad=True)
				dec=dec.to_string(sep='dms',pad=True,)
				print(ra[:5],dec[:5])

				return [ID,ra,dec]

		else: 
			
			print('A ' + self.name + ' file is missing! Abort! Abort!')
			
			return

	#--------------------------------------------------------------------------------------------------
		
	def loadcatalogs(self,debug=False):
		
		path=sys.argv[1]
		os.chdir(path)
		
		print('Currently working in \"' + os.getcwd()+ '\"')

		self.sampedro_n0=self.loadcatalog(path,'_Sampedro_cluster_members_query.csv',(2,3,5,57,58,59),'\t',1 ,[('ID','U10'),('ra','f8'),('dec','f8'),('M1','b'),('M2','b'),('M3','b')],debug=debug)
		dt = np.dtype([('ID', np.unicode_,21), ('ra', np.unicode_, 12), ('dec', np.unicode_, 11)])
		
		#self.PS=self.loadcatalog(path, '_panstarrs_search.txt',(1,2,3),'\t',2 ,dt,debug=debug)
		self.PS=self.loadcatalog(path, '_panstarrs_search.csv',(0,1,2),',',1 ,dt,debug=debug)
		self.K2=self.loadcatalog(path, '_k2_search.txt',(0,4,5,20),'\t',2 ,[('ID','i8'),('ra','U12'),('dec','U12'),('2MASS','U12')],debug=debug)
		
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
		print('Number of Everest LCs with 2MASS ID:')		
		print(len(self.K2MASS[2]))
		
		return

	#--------------------------------------------------------------------------------------------------

	def refinesampedro(self,debug=False):
	
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
				self.sampedro_n3=[[0],['0h0m0s'],['0h0m0s'],[0],[0],[0]]
			
			else:
				for _ in range(len(self.sampedro_n2)):
	
					self.sampedro_n3[_] = [i for j, i in enumerate(self.sampedro_n3[_]) if j not in dell]
		
		print('Total number of stars in Sampedro:')
		print(len(self.sampedro_n0[5]))
		print('...with one membership affirmation:')
		print(len(self.sampedro_n1[5]))
		print('...with 2 membership affirmations:')
		print(len(self.sampedro_n2[5]))
		print('...with 3 membership affirmations:')
		print(len(self.sampedro_n3[5]))

		return

	#---------------------------------------------------------------------------------------------------------------------

	def sampedro_match(self, n, dist='0d0m3s', cat='Pan-STARRS',debug=False):
			
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
		
		ra2=l[1]
		dec2=l[2]
		ra3=sampedro[n][1] 
		dec3=sampedro[n][2]

		#match Sampedro cluster members with Pan-STARRS
		try:

			match_idx, match_K2ID, match_K2Ra,match_K2Dec,len_match_idx=match(l,ra2,dec2,ra3,dec3,dist)
		
		except ValueError:

			print(cat + 'catalog list or Sampedro_n' + str(n) + ' are probably empty. No matching possible.')
			match_idx=[0]


	
		return match_idx, match_K2ID, match_K2Ra,match_K2Dec,len_match_idx

	#--------------------------------------------------------------------------------------------

	def second_order_match(self, n, distPS='0d0m3s', dist='0d0m3s',debug=False):

		'''
		
		Matching K2 and Pan-STARRS (conditional on Sampedro cluster membership with (n) shared assessments).
		
		'''

		lK2=list(self.K2)
		lPS=list(self.PS)
		#THE PROBLEM MAY BE HERE
		match_idx_PS, match_K2ID, match_K2Ra,match_K2Dec,length=self.sampedro_match(n, dist=distPS, cat='Pan-STARRS',debug=debug)

		idx_K2=list(range(len(lK2[0])))
		
		for _ in range(len(self.K2)):
			lK2[_]= [i for j, i in enumerate(lK2[_]) if j in idx_K2]

		for _ in range(len(self.PS)):
			lPS[_]= [i for j, i in enumerate(lPS[_]) if j in match_idx_PS]
		print(len(lPS[0]))
		print('From the ' + self.name + ' members as assessed by Sampedro '+str(n)+' ' + str(len(lPS[0])) + ' have an Pan-STARRS ID.')
		if len(lPS[1])==0 or len(lK2[1])==0:
			print('The arrays are empty for Sampedro_' + str(n)+'.')
			second_match_idx=[]
		else:
			ra2=lK2[1]
			dec2=lK2[2]
			
			ra3=lPS[1] 
			dec3=lPS[2]

			second_match_idx, match_K2_ID, match_K2_Ra, match_K2_Dec,len_second_match_idx =match(lK2,ra3,dec3,ra2,dec2,dist)

		return second_match_idx, match_K2_ID, match_K2_Ra, match_K2_Dec, len_second_match_idx

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------

#General purpose funcs:
def match(cat,ra1,dec1,ra2,dec2,dist='0d0m1s'):#ra2, dec2 correspond to cat
    '''
    Finds the entries in a catalog cat with corresponsing ra2 and dec2 that match with a list or coordinates (ra1,dec1)  
    '''
    ra1,ra2,dec1,dec2=Angle(ra1),Angle(ra2),Angle(dec1),Angle(dec2)
    c = ICRS(ra=ra2, dec=dec2)
    catalog = ICRS(ra=ra1, dec=dec1)
    idx, d2d, d3d = match_coordinates_sky(c, catalog,nthneighbor=1) #idx are indices into catalog that are the closest objects to each of the coordinates in c
    second_match_idx=are_within_bounds(idx, d2d,'0d0m0s', dist)
    match_ID=[i for j, i in enumerate(cat[0]) if j in second_match_idx]
    match_Ra=[i for j, i in enumerate(cat[1]) if j in second_match_idx]
    match_Dec=[i for j, i in enumerate(cat[2]) if j in second_match_idx]
			
    return second_match_idx,match_ID,match_Ra, match_Dec, str(len(second_match_idx))

def mocktest():
	
	x=OpenCluster('mockup','mockup', radius=10, age=1)
	x.loadcatalogs()
	x.refinesampedro()
	
	try:
		l, k2,lengtha=x.sampedro_match(1,dist='0d0m5s',cat='Pan-STARRS')
		l, k2,lengthb=x.sampedro_match(1,dist='0d0m3s',cat='Pan-STARRS')
		l, k2,lengthc=x.sampedro_match(1,dist='0d0m5s',cat='K2MASS')
		l, k2,lengthd=x.sampedro_match(1,dist='0d0m3s',cat='K2MASS')
		l, lengthe=x.second_order_match(1,distPS='0d0m5s',dist='0d0m5s')
		l, lengthf=x.second_order_match(1,distPS='0d0m5s',dist='0d0m3s')
		l, lengthg=x.second_order_match(1,distPS='0d0m3s',dist='0d0m5s')
		l, lengthh=x.second_order_match(1,distPS='0d0m3s',dist='0d0m3s')
	
		a=[len(x.sampedro_n0[0]),len(x.sampedro_n1[0]),len(x.PS[0]),len(x.K2[0]), lengtha,lengthb,lengthc,lengthd,lengthe,lengthf,lengthg,lengthh]
		a=[int(i) for i in a]
		assert a==[40,40,40,20,10,10,2,2,5,5,5,5]
		print('Test success. Green light for matching.')
	except AssertionError:
		print('Your results may not be valid. Test failed. Please check.')		
	return

def find(name, path,debug=False):
    
    '''
    Return the path to a file with name if it exists within the folder in path.
    
    Args:
    name: name of the file one is looking for
    path: path where on suspects to find the file
    debug: generates extra output to the command prompt for debugging
    
    Returns:
    path to file
    '''
    
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

def are_within_bounds(idx,d2d, min_angle, max_angle,debug=False):
	
	'''
	Check if the objects with indices idx which match a catalog object within distance d2d lie within [min_angle, max-angle].
	
	Args:
	idx: list of indices into said catalog for identification
	d2d: list of corresponding distances
	min_angle: minimum distance in dms format
	max_angle: maximum distance in dms format
	debug: generates extra output to the command prompt for debugging
	
	Returns:
	n: indices extracted from idx that fulfil the condition for angular distance
	'''
	
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
		if len(m)>1:
			for k in range(len(m)-1):
				if dpop[m[k]].is_within_bounds(dpop[m[k+1]]): 
					minidx=lpop[m[k]]
				else: 
					minidx=lpop[m[k+1]]
		else: minidx=lpop[m[0]]
		n.append(minidx)
		lpop=[i for j, i in enumerate(lpop) if j not in m]
		dpop=[i for j, i in enumerate(dpop) if j not in m]
	return n

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


def wrap_cross(inputs,debug=False):
	for item in inputs:
		
		x=OpenCluster(item[0],item[1], radius=item[2], age=item[3])
		x.loadcatalogs(debug=debug)
		print('\nMatching catalogs for ' + x.name + ':\n')
		x.refinesampedro(debug=debug)
		id2mass=open('share/2mass_match/'+item[1]+'_2mass_IDs.txt', 'w')
		idps=open('share/panstarrs_match/'+item[1]+'_panstarrs_IDs.txt', 'w')
		idcross=open('share/cross_match/'+item[1]+'_cross_IDs.txt', 'w')
		#out.write(str(len(x.sampedro_n0[0]))+'\n'+str(len(x.sampedro_n1[0]))+'\n'+str(len(x.sampedro_n2[0]))+'\n'+str(len(x.sampedro_n3[0]))+'\n'+str(len(x.PS[0]))+'\n'+str(len(x.K2[0]))+'\n')
		lps_m, ps_m, psra_m,psdec_m,len1=x.second_order_match(1,distPS='0d0m3s',dist='0d0m3s',debug=debug)
		lk2_m, k2_m,k2ra_m,k2dec_m, len2=x.sampedro_match(1,dist='0d0m3s',cat='K2MASS',debug=debug)
		print(ps_m,k2_m)
		cross=set(ps_m) & set(k2_m)
		for x in k2_m:
			id2mass.write('EPIC '+str(x)+'\n')
			
		for x in ps_m:
			idps.write('EPIC '+str(x)+'\n')
			
		for x in list(cross):
			idcross.write('EPIC '+str(x)+'\n')
		idps.close()
		id2mass.close()
		idcross.close()
	return
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

def wrap_basic(inputs,debug=False):
	for item in inputs:
		
		x=OpenCluster(item[0],item[1], radius=item[2], age=item[3])
		x.loadcatalogs(debug=debug)
		print('\nMatching catalogs for ' + x.name + ':\n')
		x.refinesampedro(debug=debug)
		out=open('share/'+item[1]+'_IDs.txt', 'w')
		out.write(str(len(x.sampedro_n0[0]))+'\n'+str(len(x.sampedro_n1[0]))+'\n'+str(len(x.sampedro_n2[0]))+'\n'+str(len(x.sampedro_n3[0]))+'\n'+str(len(x.PS[0]))+'\n'+str(len(x.K2[0]))+'\n')
		

		#all line numbers +2 in parameters sheet
		for i in range(1,4):
			##10
			l, k2, ra,dec,length=x.sampedro_match(i,dist='0d0m5s',cat='Pan-STARRS',debug=debug)
			out.write(length+'\n')
			##11
			l, k2, ra,dec,length=x.sampedro_match(i,dist='0d0m3s',cat='Pan-STARRS',debug=debug)
			out.write(length+'\n')
		for i in range(1,4): 
			##16,18,20
			l, k2, ra,dec,length=x.sampedro_match(i,dist='0d0m5s',cat='K2MASS',debug=debug)
			out.write(length+'\n')
			##17,19,21
			l, k2, ra,dec,length=x.sampedro_match(i,dist='0d0m3s',cat='K2MASS',debug=debug)
			out.write(length+'\n')
		
			# if i==1:		
			#     for item in ps:
			#         ids.write(str(item)+'\n')
		for i in range(1,4):
			##22,26,30
			l, ps, ra,dec,length=x.second_order_match(i,distPS='0d0m5s',dist='0d0m5s',debug=debug)
			out.write(length+'\n')
			##23,27,31
			l, ps, ra,dec,length=x.second_order_match(i,distPS='0d0m5s',dist='0d0m3s',debug=debug)
			out.write(length+'\n')
			##24,28,32
			l, ps, ra,dec,length=x.second_order_match(i,distPS='0d0m3s',dist='0d0m5s',debug=debug)
			out.write(length+'\n')
			##25,29,33
			l, ps, ra,dec,length=x.second_order_match(i,distPS='0d0m3s',dist='0d0m3s',debug=debug)
			out.write(length+'\n')
			
			#out.write(length+'\n')
		out.close()
	return
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#Wrapping it:


inputs=[]
#inputs.append(['M67','M67', 15, 4.0])
#inputs.append(['Ruprecht 147','Ruprecht_147', 30, 2.5])
inputs.append(['M44','M44', 47, 0.73])
wrap_cross(inputs,debug=True)


