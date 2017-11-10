"""Module slicingAndSubsetting

Contains functions to extract subsets of initial data (flattening the data or 
retaining the dimensions for further operations), or execute some operations
on them while retaining the other dimensions.
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import dask.array as da
from netCDF4 import Dataset
from scipy.interpolate import interp1d

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from daskOptions import *
    
#---- Parameters ----#

#---- Functions ----# 

## Extracts a subset of data points for further statistical analysis.
def sample1d(values,subset=None):

	"""Arguments:
		- values: a 3D (time,lat,lon) numpy array
		- subset: a 2D (lat,lon) or 3D (time,lat,lon) numpy array of bools
	Returns:
		- 1D numpy.array with valid data"""

	# Reduce values to subset
	if subset is None:	# Then use all points
		n_pts = values.size
		return values.flatten()
	# Otherwise, take slices at each time step
	n_pts = subset.sum()
	sshape = subset.shape
	if len(sshape) == 2:
		ntime = values.shape[0]
		ind_sub = np.vstack([[subset]]*ntime)
	elif len(sshape) == 3:
		ind_sub = subset

	return values[ind_sub]

## Computes mean of values over selected points
def computeTimeHorizontalMean(values,ind,is_3D,da_compute=da_compute_default):

	"""Arguments:
		- values (dask.array or numpy.ndarray)
		- ind (numpy or dask). It must have dimensions consistent with 'values':
		if they have same shape, their dimensions must be the same. Otherwise,
		the only allowed missing dimensions is the vertical, and if so all other
		dimensions must match.
	Returns:
		- Computes the mean over selected points. Retains the vertical dimension
		if it exists.
	Use:
		- Computes the mean over points corresponding to a given percentile bin
		of some distribution."""

	cn = getArrayType(values)
	vshape = values.shape
	ndims = len(vshape)
	ishape = ind.shape

	# Extend dimensions of indices
	if is_3D and ishape != vshape:
		# Then check that all dimensions match except the vertical one (which 
		# must be in second position)
		if vshape[0] != ishape[0] or vshape[2:] != ishape[1:]:
			print("Error in computeMeanAtTimeLatLonIndices:")
			print("    Sizes don't match along time-horizontal dimensions.")
			return
		# Extend indices along a new vertical axis
		nlev = vshape[1]
		ind = cn.repeat(ind[:,np.newaxis,...],nlev,1)

	# Compute time-horizontal mean
	v_nan = np.copy(values)
	v_nan[ind != 1] = np.nan
	if is_3D:
		nlev = vshape[1]
		dimorder = [1,0]+list(range(2,ndims))
		v_mean = cn.nanmean(cn.transpose(v_nan,dimorder).reshape(nlev,-1),axis=1)
	else:
		v_mean = cn.nanmean(v_nan)

	# Compile operations within function if necessary
	if ind.__class__ == da.core.Array and da_compute:
		v_mean = v_mean.compute()

	return v_mean


## Compute time-horizontal mean over all points
def computeMean(values,subset=None):

	"""-- Yet to implement. See functions/statistics.py in SPCAM Extremes.--"""


## Reduce domain to prescribed area
def reduceDomain(values,subsetname,examplefile,varid):

	"""Extracts a subset of points while retaining the object shape. Typically,
	take a slice in a given dimension.
		- values can be numpy.ndarray or dask.array. 
		- examplefile is used to get the size and index of the dimension of
		interest.
		-subsetname indicates which points to extract (e.g. 'tropics')"""

	cn = getArrayType(values)

	if subsetname is not 'global':
		fh = Dataset(examplefile,'r')
		dims = fh.variables[varid].dimensions
		ndims = len(dims)
		fh.close()
	if subsetname == 'tropics':
		indlat = np.argmax(np.array(dims) == 'lat')
		nlat = values.shape[indlat]
		lat_slice = list(range(nlat//3,2*(nlat//3)))
		values = cn.take(values,lat_slice,axis=indlat)

	return values

## Get isobaric surfaces
def isobaricSurface(var,pres,p_ref=500,levdim=1):
	
	"""From a 4D variables, extract a 3D (time-lat-lon) surface at constant 
	pressure p_ref (in hPa).
	Assumes values of p are decreasing along levdim."""

	cn = getArrayType(var)

	p_ref *= 100	# in Pa
	# Put vertical dimension in last position
	dimorder = list(range(len(pres.shape)))
	newdimorder = dimorder[:levdim]+dimorder[levdim+1:]+[dimorder[levdim]]
	pres_perm = cn.transpose(pres,newdimorder)
	# Define new shape
	var_ravel = cn.ravel(cn.transpose(var,newdimorder))
	newshape = list(pres.shape)
	del newshape[levdim]
	# Define new chunks
	newchunks = None
	if cn == da:
		newchunks = pres.chunks[:levdim]+pres.chunks[levdim+1:]
	# Ravel pressure coordinate
	pres_ravel = cn.ravel(pres_perm)
	# Extract indices below (i.e. larger than) p_ref
	i_LT = (pres_ravel >= p_ref).astype(int)
	# Ravelled indices just below and just above reference pressure
	i_diff = cn.diff(i_LT)
	i_diff[i_diff == -1] = 0
	i_diff2 = np.take(i_LT,range(i_LT.size-1),axis=-1) - np.take(i_LT,range(1,i_LT.size),axis=-1)
	i_diff2[i_diff2 == -1] = 0
	i_below = cn.hstack((np.zeros((1,),dtype=int),i_diff))
	i_above = cn.hstack((i_diff,np.zeros((1,),dtype=int)))
	i_below,i_above = i_below.astype(bool,copy=False),i_above.astype(bool,copy=False)
	# Compute fraction coefficient to interpolate values
	i_valid = np.mean((pres_perm >= p_ref),axis=-1,dtype=bool).flatten()
	
	def getVals(vals_ravel,i_valid,stencil,newshape,newchunks=None):
		values = np.nan*np.ones(i_valid.shape)
		values[i_valid] = vals_ravel[stencil]
		if cn == da:
			return da.from_array(np.reshape(values,newshape),chunks=newchunks)
		elif cn == np:
			return np.reshape(values,newshape)

	pres_below = getVals(pres_ravel,i_valid,i_below,newshape,newchunks)
	pres_above = getVals(pres_ravel,i_valid,i_above,newshape,newchunks)

	# if cn == da:
	# 	pres_below = da.compress(i_below.compute(),pres_ravel)
	# 	pres_above = da.compress(i_above.compute(),pres_ravel)
	# elif cn == np:
	# 	pres_below = getVals(pres_ravel,i_valid,i_below,newshape)
	# 	pres_above = getVals(pres_ravel,i_valid,i_above,newshape)


	f = (p_ref-pres_above)/(pres_below-pres_above)
	# Interpolate variable onto p_ref surface
	
	# if cn == da:
	# 	var_pref = cn.reshape(f*da.compress(i_below.compute(),var_ravel) + 
	# 						  (1-f)*da.compress(i_above.compute(),var_ravel),newshape)
	# elif cn == np:
	# 	# var_pref = cn.reshape(f*var_ravel[i_below] + (1-f)*var_ravel[i_above],newshape)
	
	var_pref = f*getVals(var_ravel,i_valid,i_below,newshape,newchunks) +\
				   (1-f)*getVals(var_ravel,i_valid,i_above,newshape,newchunks)

	return var_pref

## Other version to interpolate variable on pressure level
def varAtPressureLevelInterp1D(var,pres3D,p_ref):

	"""Only works with numpy arrays, not dask.array's.
	After implementing the same function using scipy.interpolate.griddata,
	this is substantially faster."""

	n = pres3D[0].size
	pshape = pres3D[0].shape
	var_ref = np.array([np.nan]*n).reshape(pshape)

	for (ilat,ilon) in np.ndindex(*pshape):
		try:
			var_ref[ilat,ilon] = interp1d(pres3D[:,ilat,ilon],var[:,ilat,ilon])(p_ref)
		except ValueError:
			var_ref[ilat,ilon] = np.nan

	return var_ref

## Interpolates var at a given pressure level
def varAtPressureLevel(var,pres3D,p_ref,timedim=0,levdim=1):

	"""This is a wrapper around varAtPressureLevelInterp1D to work with numpy
	or dask arrays. Careful, varAtPressureLevelInterp1D assumes a specific order
	timedim,levdim,latdim,londim."""

	cn = getArrayType(var)
	vshape = var.shape

	out_list = []
	for itime in range(vshape[timedim]):
		v = cn.squeeze(cn.take(var,[itime],axis=timedim),axis=timedim)
		p = cn.squeeze(cn.take(pres3D,[itime],axis=timedim),axis=timedim)
		if cn == da:
			v,p = v.compute(), p.compute()
		v_out = varAtPressureLevelInterp1D(v,p,p_ref)
		out_list.append(v_out[np.newaxis,...])
	var_out = np.vstack(out_list)

	if cn == da:
		newchunks = list(var.chunks)
		del newchunks[levdim]
		var_out = da.from_array(var_out,chunks=tuple(newchunks))

	return var_out



