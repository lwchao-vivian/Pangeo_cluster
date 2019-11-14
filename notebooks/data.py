import numpy as np
import scipy.interpolate
from netCDF4 import Dataset
import pdb

class Data3d:
	"""basic data element that contains a 3D (plevs/lat/lon) data field"""

	def __init__(self,array3d=[],lon=[],lat=[],plevs=[],time=[],minv=-9e9):
		"""Data3d(array[time,plevs,lat,lon], lon, lat, plevs, time,minv): 3D data field object time is optional.
		Values less than minv are masked."""
		
		# if called w/ no parameters
		if np.size(array3d) == 0:
			self.data=[];self.lon=[]
			self.lat=[];self.time=[]
			self.plevs=[]
			return
		
		if len(np.shape(array3d)) == 2:
			array3d=np.reshape(array3d,(1,1,len(lat),len(lon)))
		if len(np.shape(array3d)) == 3:
			if len(time) > 0 and len(plevs) == 0:
				array3d=np.reshape(array3d,(len(time),1,len(lat),len(lon)))
			if len(time) == 0 and len(plevs) > 0:
				array3d=np.reshape(array3d,(1,len(plevs),len(lat),len(lon)))
			if len(time) == 0 and len(plevs) == 0:
				array3d=np.reshape(array3d,(1,1,len(lat),len(lon)))
		if len(np.shape(array3d)) != 4:
			raise TypeError("requires 3D or 4D input data")
		s1=np.shape(array3d)
		if len(plevs) > 0 and len(plevs) != s1[1]:
			raise TypeError("pressure mis-match")
		if len(lat) != s1[2]:
			raise TypeError("longitude mis-match")
		if len(lon) != s1[3]:
			raise TypeError("latitude mis-match")

		# reverse latitudes if descending
		if lat[0] > lat[1]:
			s1=np.shape(array3d)
			for num in range(s1[0]):
				for num2 in range(s1[1]):
					x1=np.ma.copy(array3d[num,num2,])
					x1=np.flipud(x1)
					array3d[num,num2,]=x1
			lat=np.flipud(lat)

		# convert -180-180 to 0-360 longitudes
		if min(lon) < -50:
			x1=len(lon)/2
			lon=np.concatenate( (lon[x1:],lon[0:x1]+360) )
			array3d=np.ma.concatenate( (array3d[:,:,:,x1:],array3d[:,:,:,0:x1]), axis=3 )

		array3d=np.ma.array(array3d) # mask bad values
		np.ma.masked_where(array3d < minv,array3d,copy=False)
		np.ma.masked_where(np.isnan(array3d),array3d,copy=False)
		
		self.data=array3d
		self.lat=np.array(lat)
		self.lon=np.array(lon)
		self.plevs=np.array(plevs)
		self.time=np.array(time)
		return

	def copy(self):
		"""produce a copy of the Data3d structure"""
		newdata=Data3d(np.ma.copy(self.data),np.copy(self.lon),np.copy(self.lat),np.copy(self.plevs),np.copy(self.time))
		return newdata

	def append(self,newdata):
		"""appends time slice(s) to the 3D time series"""

		if np.size(self.data) == 0:
			self.data=newdata.data
			self.lon=newdata.lon
			self.lat=newdata.lat
			self.plevs=newdata.plevs
			self.time=newdata.time
			return

		z1=np.append(self.data,newdata.data,axis=0)
		self.data=z1
		z2=np.append(self.time,newdata.time)
		self.time=z2

		return 
 
	def zonal_average(self):
		"""zonal average of data"""
		# NOT WORKING
		xxx
		return np.ma.average(self.data,axis=3)

	def global_averageOLD(self,latrange=(-90,90),minv=-9e9):
		"""x1.global_average( latrange=(minlat, maxlat), minv=minv )
		average over latitude range"""
		
		lat=self.lat
		ind=(lat >= min(latrange)) & (lat <= max(latrange))
		lat=lat[ind]
		d1=self.data[:,:,ind,]
		
		coslat = np.cos(lat*3.14159/180.)
		c1=coslat.reshape(len(coslat),1)
		wgt=np.zeros( (len(lat),len(self.lon)) )+c1

		ga=np.ma.zeros( (max(1,len(self.time)),max(1,len(self.plevs))) )
		for num in range(max(1,len(self.time))):
			for pnum in range(max(1,len(self.plevs))):
				d2=d1[num,pnum,]
				ind2=np.ma.where(d2 <= minv)
				wgt1=wgt;wgt1[ind2]=0 
				# ga[num,pnum]=np.ma.average(self.data[num,pnum,ind,:],weights=wgt)
				ga[num,pnum]=np.ma.average(d2,weights=wgt1)
				
		if min(np.shape(ga)) == 1:
			ga=np.copy(ga.flatten())
			
		return ga

	def global_average(self,latrange=(-90,90)):
		"""x1.global_average( latrange=(minlat, maxlat) )
		average over latitude range"""
		
		coslat = np.cos(self.lat*3.14159/180.)
		ind=(self.lat >= min(latrange)) & (self.lat <= max(latrange))

		x1=np.ma.average(self.data,axis=3)*coslat.reshape(1,1,len(coslat))
			
		return np.squeeze(np.ma.sum(x1[:,:,ind],axis=2)/np.ma.sum(coslat[ind]))

	def anomaly(self,timerange=(-9e9,9e9),minv=-9e9):
		"""x1.anomaly(timerange=(mintime,maxtime))
		calculate anomaly relative to a particular time span """

		d1=np.ma.copy(self.data)
		ind=np.where((self.time >= min(timerange)) & (self.time <= max(timerange)))[0]
		
		indx=np.arange(len(self.time)) % 12
		for ii in range(12):
			xx=np.where(ii == indx)[0] # index into data
			yy=np.where(ii == indx[ind])[0] # index into reference period 
			d1[xx,]=d1[xx,]-np.ma.average(d1[ind[yy],],axis=0)
			
		ga=Data3d(d1,self.lon,self.lat,self.plevs,time=self.time)
		
		return ga

	def interpolate(self,newlon,newlat,newplevs=[]):
		"""interpolate(newlon, newlat, newplevs)"""
	   
		newdata=np.zeros( (max(1,len(self.time)),max(1,len(self.plevs)),len(newlat),len(newlon)) )
		if len(newplevs) == 0: newplevs = self.plevs
		
		for num in range(max(1,len(self.time))):
			for num2 in range(max(1,len(self.plevs))):
				spl = scipy.interpolate.RectBivariateSpline(self.lat,self.lon,
					self.data[num,num2,]) #,kx=1,ky=1) 
				n1 = spl(newlat,newlon)
				newdata[num,num2,]=n1

		if len(newplevs) + len(self.plevs) == 0:
			return Data3d(newdata,newlon,newlat,self.plevs,self.time)

		if len(newplevs) == len(self.plevs):
			if max(newplevs-self.plevs) == 0:
				return Data3d(newdata,newlon,newlat,self.plevs,self.time)

		newdata2=np.zeros( (len(self.time),len(newplevs),len(newlat),len(newlon)) )
		newp=map(np.log,newplevs)
		oldp=map(np.log,self.plevs)
		if np.all(np.diff(oldp) > 0): # test if plevs are increasing
			for num in range(len(newlon)):
				for num2 in range(len(newlat)):
					for num3 in range(len(self.time)):
						newdata2[num3,:,num2,num]=np.interp(newp,oldp,newdata[num3,:,num2,num])
		else:
			for num in range(len(newlon)): # loop for decreasing plevs
				for num2 in range(len(newlat)):
					for num3 in range(len(self.time)):
						newdata2[num3,:,num2,num]=np.interp(newp,np.flipud(oldp),np.flipud(newdata[num3,:,num2,num]))
	
		return Data3d(newdata2,newlon,newlat,newplevs,self.time)
	def timeclip(self,r):
		"""r is a 2-element list; this routine clips the data set so that the new time falls into the range r
		(note: r could actualy be a bigger array; in that case, it uses the min and max of r to do the limit)"""

		ind=np.where((np.min(r) <= self.time) & (self.time <= np.max(r)))[0]
		# if copy == False:
		#	  self.time=self.time[ind]
		#	  self.data=self.data[ind,]
		#	  return

		new3d=self.copy()
		new3d.time=self.time[ind]
		new3d.data=self.data[ind,]	  
		return new3d

class TimeSeries:
	"""basic data element that contains a 2D string of data (time, pressure) and a time variable"""

	def __init__(self,array3d=[],time=[],plevs=[],minv=-9e9):
		"""TimeSeries(array[time], time, plevs, minv): time series object
		time is optional, filled in increasing integers if omitted.	 
		Values less than minv are masked."""

		# if called w/ no parameters
		if np.size(array3d) == 0:
			self.data=[];self.time=[];self.plevs=[]
			return

		if len(np.shape(array3d)) > 2:
			raise TypeError("requires 1D or 2D input data")
		if len(time) == 0: time=np.arange(np.shape(array3d)[0])

		if len(plevs) == 0:
			plevs=[0]
			if len(np.shape(array3d)) > 1: plevs=np.arange(np.shape(array3d)[1])
			
		if np.shape(array3d)[0] != len(time):
			raise TypeError("time array mismatch")
			
		array3d=np.ma.array(array3d) # mask bad values
		array3d=np.ma.masked_where(array3d < minv,array3d)
		array3d=np.ma.masked_where(np.isnan(array3d),array3d)

		self.data=array3d
		self.time=np.array(time)
		self.plevs=np.array(plevs)
		return

	def copy(self):
		"""produce a copy of the Data3d structure"""
		newdata=TimeSeries(np.ma.copy(self.data),np.ma.copy(self.time))
		return newdata

	def append(self,newdata):
		"""appends time slice(s) to the 3D time series"""

		if np.size(self.data) == 0:
			self.data=newdata.data
			self.time=newdata.time
			return

		if len(np.shape(self.data)) > 1:
			if np.shape(self.data)[1] != np.shape(self.data):
				raise TypeError("shape error")

		z1=np.append(self.data,newdata.data,axis=0)
		self.data=z1
		z2=np.append(self.time,newdata.time)
		self.time=z2

		return 

	def anomaly(self,timerange=(-9e9,9e9),minv=-9e9):
		"""x1.anomaly(timerange=(mintime,maxtime))
		calculate anomaly relative to a particular time span """

		d1=np.ma.copy(self.data)
		time=np.copy(self.time)
		ind=np.where((time >= min(timerange)) & (time <= max(timerange)))[0]

		indx=np.arange(len(time)) % 12
		if len(np.shape(self.data)) > 1: # lp2 is the number of plevs
			lp2=np.shape(self.data)[1]
		else:
			lp2=1
		for ii in range(12):
			for jj in range(lp2):
				xx=np.where(ii == indx)[0] # index into data
				yy=np.where(ii == indx[ind])[0] # index into reference period 
				if lp2 == 1:
					d1[xx]=d1[xx]-np.ma.average(d1[ind[yy]])
					# d1[xx]=(d1[xx]-np.ma.average(d1[ind[yy]]))/np.std(d1[ind[yy]])
				else:
					d1[xx,jj]=d1[xx,jj]-np.ma.average(d1[ind[yy],jj])
					# d1[xx,jj]=(d1[xx,jj]-np.ma.average(d1[ind[yy],jj]))/np.std(d1[ind[yy],jj]) # removes StDev and average

		ga=TimeSeries(d1,time,self.plevs)

		return ga

	def interpolate(self,newtime):
		"""interpolate(newtime)"""
		
		if len(np.shape(self.data)) > 1: # handle multiple levels
			newd=[]
			for ii in range(np.shape(self.data)[1]):
				newd.append(np.interp(newtime,self.time,self.data[:,ii],left=np.nan,right=np.nan))
			newd=np.array(newd).T
		else: # if there's only one level
			newd=np.interp(newtime,self.time,self.data,left=np.nan,right=np.nan)

		return TimeSeries(newd,newtime,self.plevs)

	def lag(self,newtime):
		"""lag(lag): positive value moves the time series forward in time"""
		
		x1=self.interpolate(self.time-newtime)
		x1.time=self.time
		return x1

	def timeclip(self,r):
		"""r is a 2-element list; this routine clips the data set so that the new time falls into the range r
		(note: r could actualy be a bigger array; in that case, it uses the min and max of r to do the limit)"""
		ind=np.where((np.min(r) <= self.time) & (self.time <= np.max(r)))[0]
		new3d=TimeSeries(self.data[ind,],self.time[ind],self.plevs)
		return new3d

	# def skip(self,r):
	#	  """r is an integer; this returns a TimeSeries with every rth element"""
	#	  x1=len(self.data)/r
	#	  ind=np.arange(x1)*r
	#	  new3d=TimeSeries(self.data[ind],self.time[ind])
	#	  return new3d

#########

def readncfile(fn,dataname,height=[]):
	"""readncfile(fn, dataname, height): read nc file with file name fn, 
	dataname name data, and height named height"""
	fni=Dataset(fn)
	k1=[str(x) for x in fni.variables.keys()]
	k1=map(str.upper,k1)

	# get longitude variable
	loni=-1;lati=-1;timi=-1;hei=-1
	
	num=0
	for k1st in k1:
		# if str.find(k1st,'LON') >= 0:
		if k1st[0:3]=='LON' and k1st.find('BNDS') == -1:
			loni=num
		# if str.find(k1st,'LAT') >= 0:
		if k1st[0:3]=='LAT' and k1st.find('BNDS') == -1:
			lati=num
		if k1st == 'TIME' or k1st == 'DATE':
			timi=num
		if len(height) > 0:
			# if str.find(k1st,height.upper()) >= 0 and k1st.find('BNDS') == -1:
			if k1st[:len(height)] == height.upper() and k1st.find('BNDS') == -1:
				hei=num
		num += 1

	if loni == -1:		  
		raise TypeError("no longitude variable")
	if lati == -1:		  
		raise TypeError("no latitude variable")
   
	# read variables
	k1=fni.variables.keys()
	d1=np.ma.array(readvar(fni,dataname))
	lon=np.array(readvar(fni,k1[loni]))
	lat=np.array(readvar(fni,k1[lati]))
		
	if hei != -1:
		plevs=np.array(readvar(fni,k1[hei]))
	else:
		plevs=[]
	if timi != -1:
		tim=np.array(readvar(fni,k1[timi]))
	else:
		tim=[]
	fni.close()

	sh=np.shape(d1)
	if len(lon) != sh[-1]: # zonal avg. file
		lon=[180.];d1=d1.reshape(sh[0],sh[1],sh[2],1)
		
	# return d1,lon,lat,plevs
	return Data3d(d1,lon,lat,plevs,tim)

def readvar(fni,varname,latrange=[],prange=[]):
	"""reads variable from netcdf object in varname
	fni is a Dataset variable"""
	tsv = fni.variables[varname]
	ts=tsv[:]
	
	# try:
	#	  ts = np.array(ts,dtype='float32')*float(tsv.scale_factor)+float(tsv.add_offset)
	# except AttributeError:
	#	  x=0 # dummy statement

	try:
		if tsv.calendar == 'noleap': ts /= 365.
		if tsv.calendar == '365_day': ts /= 365.
		if tsv.calendar == '360_day': ts /= 360.
		if tsv.calendar == 'gregorian': ts /= 365.24
	except AttributeError:
		x=0 # dummy statement

	return ts

def wvsat(t,p):
	"""wvsat(t(K), p(mb)): return water vapor saturation mixing ratio (g/kg)
	over liquid for t > 273 K and over ice for t < 273 K"""
	
	if type(t) is float:
		return float(wvsat(np.array([t]),np.array([p])))
	
	Psat=t*0.

	ind1=np.ma.where(t >= 273.15)
	ind2=np.ma.where(t < 273.15)
	
# Source Hyland, R. W. and A. Wexler, Formulations for the Thermodynamic Properties of the saturated Phases of H2O from 173.15K to 473.15K, ASHRAE Trans, 89(2A), 500-519, 1983.
	Psat[ind1] = np.exp(  -0.58002206e4 / t[ind1] + 0.13914993e1 - 0.48640239e-1 * t[ind1] 
	+ 0.41764768e-4 * t[ind1]**2. - 0.14452093e-7 * t[ind1]**3. 
	+ 0.65459673e1 * np.log(t[ind1]) ) / 100.
	
# Source : Goff-Gratch, Smithsonian Meteorological Tables, 5th edition, p. 350, 1984
	ei0	   = 6.1071		  # mbar
	T0	   = 273.16		  # freezing point in K
	Psat[ind2] = 10.**(-9.09718 * (T0 / t[ind2] - 1.) - 3.56654 * np.log10(T0 / t[ind2]) + 0.876793 * (1. - t[ind2] / T0) + np.log10(ei0))
	
	Psat = Psat/p*18/29.*1e3  # convert to g/kg
	return Psat

def anomaly(data):
	"""calculate anomaly of a 1D array (assumes no missing time slots)"""
	d1=np.ma.copy(data)
	ind=np.ma.array(range(len(data))) % 12
	for num in range(12):
		d1[num == ind]=data[num == ind]-np.ma.average(data[num == ind])
	return d1

def check(var):
	"""print out statistics of variable"""
	print("type: ",type(var))
	print("shape: ",np.shape(var))
	print("min/max: ",np.ma.min(var)," ",np.ma.max(var))

	return

def global_average(idata,lat):
	"""global_average(data,lat): average over latitude range
	data([time], lat,lon) is a 2D or 3D array, and lat is a 1D array"""

	def global_average1(idata,lat):
		"""global_average(data,lat): average over latitude range
		data(lat,lon) is a 2D array, and lat is a 1D array"""

		coslat = np.cos(lat*3.14159/180.)
		c1=coslat.reshape(len(coslat),1)
		wgt=np.zeros( idata.shape )+c1
		ga=np.ma.average(idata,weights=wgt)

		return ga

	if len(np.shape(idata)) == 2: return global_average1(idata,lat)

	sh=np.shape(idata)
	ga=[global_average1(idata[ii,],lat) for ii in range(sh[0])] 
	return np.array(ga)
