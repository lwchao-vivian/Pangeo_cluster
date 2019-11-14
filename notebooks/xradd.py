## routines to use with xarray

import numpy as np
import xarray as xr

def fixdims(var):
	"""changes the dimensions to lat, lon, and level"""
	
	for ii in var.dims:
		kk=ii[:3].lower()
		
		if kk == 'lat':
			var=var.rename({ii:'lat'})
			
		if kk == 'lon':
			var=var.rename({ii:'lon'})
			
		if kk == 'lev' or kk == 'ple' or kk == 'pre':
			var=var.rename({ii:'level'})

		if kk == 'tim':
			var=var.rename({ii:'time'})

	return var

def gavg(idata):
	"""calculate global average
	e.g., x1=gavg(d1['t2m'])"""
	
	wgt1=np.cos(np.deg2rad(idata.lat))*(idata*0+1)
	ga=(wgt1*idata).sum(dim=['lat','lon'])/wgt1.sum(dim=['lat','lon'])

	return ga

def anomaly(idata):
	"""calculate anomaly"""
	
	clim=idata.groupby('time.month').mean(dim='time')
	anom=idata.groupby('time.month')-clim
	
	return anom

# average the data to seasonal or annual time frames
seasonal=lambda ivar: ivar.resample(time='Q-FEB').mean(dim='time')
annual=lambda ivar: ivar.resample(time='A').mean(dim='time')

def covmat(a,b):
	"""covmat(a,b): calculate covariance between vector a and each grid point of matrix b
	 returns a matrix with dimensions of b"""
	
	b1=b.values
	s1=np.shape(b1)
	a1=a.values
	b1=b1.reshape(s1[0],s1[1]*s1[2])
		
	a1=a1-np.average(a1);b1=b1-np.average(b1,axis=0)
	
	c1=np.matrix(a)*np.matrix(b1)
	c1=c1.reshape(s1[1],s1[2])/s1[0]
	c1=xr.DataArray(c1,dims=['lat','lon'],coords={'lat': b.lat, 'lon': b.lon})
	
	return c1

def detrend1d(ovar):
	"""linear detrend of 1D vector; returns anomaly"""
	
	ovar1=anomaly(ovar.squeeze())

	fitx=stats.linregress(np.arange(len(ovar1)),ovar1)
	ovar1 = ovar1 - np.arange(len(ovar1))*fitx.slope

	ovar1=anomaly(ovar1)
	
	return ovar1 # return numpy array

def detrend(ovar):
	"""linear detrend of 2D data set (lat x lon x time)
	returns anomaly"""
	
	ovar1=anomaly(ovar)
	
	t1=c1=xr.DataArray(np.arange(len(ovar1.time)),dims='time',coords={'time': ovar1.time})
	slope=covmat(t1,ovar1)/np.std(t1)**2
	
	ovar1 -= slope*t1 # remove linear trend
	ovar2=anomaly(ovar1)
	
	return ovar2 