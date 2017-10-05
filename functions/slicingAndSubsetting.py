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

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from defineDaskParameters import *
    
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
def computeTimeHorizontalMean(values,ind,is_3D,da_compute=False):

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

	if values.__class__ == np.ndarray:
		cn = np # stands for 'class name'
	elif values.__class__ == da.core.Array:
		cn = da
	else:
		print("Unvalid data type:", type(sample))
		return

	vshape = values.shape
	ndims = len(vshape)
	ishape = ind.shape

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

