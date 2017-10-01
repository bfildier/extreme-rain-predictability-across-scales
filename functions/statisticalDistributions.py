"""Module statisticalDistributions

Contains functions to compute probability densities with various transformations
on the x-axis (logarithmic, linear and inverse-logarithmic).
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import dask.array as da
from scipy.interpolate import CubicSpline

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from defineDaskParameters import *
    
#---- Paramters ----#
## Default number of bins per logarithmic decade.
nbpd = 10
## Default minimum number of data points per bin.
nppb = 4	

#---- Functions ----# 

## Percentile ranks regularly spaced on an inverse-logarithmic axis (zoom on 
## largest percentiles of the distribution).
def getInvLogRanks(n_pts,n_pts_per_bin=nppb,n_bins_per_decade=nbpd):

	"""Arguments:
		- n_pts: sample size
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
		- n_bins_per_decade: number of ranks/bins per logarithmic decade
	Returns:
		- 1D numpy.array"""

	# k indexes bins
	n_decades = log10(n_pts/n_pts_per_bin) 	# Maximum number of decades
	dk = 1/n_bins_per_decade			 
	k_max = int(n_decades/dk)*dk 			# Maximum bin index
	max_rank = (1.-10**(-k_max))*100
	scale_invlog = np.arange(0,k_max+dk,dk)
	ranks_invlog = np.subtract(np.ones(scale_invlog.size),np.power(10,-scale_invlog))*100

	return max_rank, ranks_invlog

## Extracts a subset of data points for further analysis.
def getSample1D(values,subset=None):

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

## Compute percentiles of the distribution and histogram bins from percentile ranks.
def computePercentilesAndBinsFromRanks(sample,ranks):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- ranks: 1D array of floats between 0 and 1
	Returns:
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""

	if isinstance(sample,np.ndarray):
		centers = np.percentile(sample,ranks)
	elif isinstance(sample,da.core.Array):
		centers = da.percentile(sample,ranks).compute()
	breaks = np.convolve(centers,[0.5,0.5],mode='full')
	breaks[0] = CubicSpline([1,2,3],breaks[1:4],extrapolate=True)(0)
	breaks[-1] = CubicSpline([-4,-3,-2],breaks[-4:-1],extrapolate=True)(-1)
	
	return centers, breaks

## Defines percentiles and histogram bins on inverse-logarithmic ranks.
def definePercentilesOnInvLogQ(sample,n_pts_per_bin=nppb,n_bins_per_decade=nbpd):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
		- n_bins_per_decade: number of ranks/bins per logarithmic decade
	Returns:
		- max_ranks (maximum percentile rank)
		- ranks
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""
	
	n_pts = sample.size
	max_rank, ranks_invlog = getInvLogRanks(n_pts,n_pts_per_bin,n_bins_per_decade)

	centers_invlog, breaks_invlog = \
		computePercentilesAndBinsFromRanks(sample,ranks_invlog)

	return max_rank, ranks_invlog, centers_invlog, breaks_invlog

## Returns minimum and maximum values.
def findExtremumValues(sample,vmin=None,vmax=None):

	"""Argument: 1D numpy or dask array of data values."""

	# Find minimum value
	if vmin is None:
		if isinstance(sample,np.ndarray):
			vmin = np.nanmin(sample[sample > 0])
		elif isinstance(sample,da.core.Array):
			vmin = da.nanmin(sample[sample > 0])
	# Find maximum value
	if vmax is None:
		if isinstance(sample,np.ndarray):
			vmax = np.nanmax(sample)
		elif isinstance(sample,da.core.Array):
			vmax = da.nanmax(sample)

	return vmin,vmax

## Define logarithmic bin centers and edges from sample values.
def defineLogBins(sample,n_bins_per_decade=nbpd,vmin=None,vmax=None):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- n_bins_per_decade: number of ranks/bins per logarithmic decade
		- vmin and vmax: extremum values
	Returns:
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""

	vmin,vmax = findExtremumValues(sample,vmin,vmax)
	kmin = floor(log10(vmin))
	kmax = ceil(log10(vmax))
	dk = 1/n_bins_per_decade
	exps = np.arange(kmin,kmax+dk,dk)
	breaks = np.power(10.,exps)
	centers = np.convolve(breaks,[0.5,0.5],mode='valid')

	return centers, breaks

## Define linear bin centers and edges from sample values.
def defineLinearBins(sample,n_bins=50,vmin=None,vmax=None):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- vmin and vmax: extremum values
	Returns:
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""

	vmin,vmax = findExtremumValues(sample,vmin,vmax)
	breaks = np.linspace(vmin,vmax,n_bins)
	centers = np.convolve(breaks,[0.5,0.5],mode='valid')

	return centers, breaks


## Returns percentile ranks corresponding to percentile values ('centers').
def computePercentileRanksFromBins(sample,centers):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- centers: corresponding percentiles, or bin centers
	Returns:
		- ranks: 1D numpy.ndarray"""
	
	n_pts = sample.size
	ranks = np.array(list(map(lambda x:(sample < x).sum()/n_pts, centers)))

	return ranks

## Derive bins centers (percentiles) and edges, percentile ranks and corresponding 
## probability densities
def getUnivariateStatistics(sample,n_pts_per_bin=nppb,mode='linear'):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
		- mode: 'linear', 'log', 'invlogQ' 
	Returns:
		- """






