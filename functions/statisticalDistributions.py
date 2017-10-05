"""Module statisticalDistributions

Contains functions to compute probability densities with various transformations
on the x-axis (logarithmic, linear and inverse-logarithmic).
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import dask.array as da
from scipy.interpolate import lagrange

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))
    
#---- Parameters ----#
## Default number of bins per logarithmic decade.
nbpd = 10
## Default minimum number of data points per bin.
nppb = 4	
## Default number of bins used for linear statistics
nlb = 50

#---- Functions ----# 

## Percentile ranks regularly spaced on an inverse-logarithmic axis (zoom on 
## largest percentiles of the distribution).
def getInvLogRanks(n_pts,n_pts_per_bin=nppb,n_bins_per_decade=nbpd):

	"""Arguments:
		- n_pts: sample size
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
		- n_bins_per_decade: number of ranks/bins per logarithmic decade
	Returns:
		- float, 1D numpy.array"""

	# k indexes bins
	n_decades = log10(n_pts/n_pts_per_bin) 	# Maximum number of decades
	dk = 1/n_bins_per_decade			 
	k_max = int(n_decades/dk)*dk 			# Maximum bin index
	scale_invlog = np.arange(0,k_max+dk,dk)
	ranks_invlog = np.subtract(np.ones(scale_invlog.size),np.power(10,-scale_invlog))*100

	return ranks_invlog

## Compute percentiles of the distribution and histogram bins from percentile ranks.
def computePercentilesAndBinsFromRanks(sample,ranks):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- ranks: 1D array of floats between 0 and 1
	Returns:
		- ranks, cropped by one at beginning and end
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""

	if isinstance(sample,np.ndarray):
		centers = np.percentile(sample,ranks)
	elif isinstance(sample,da.core.Array):
		centers = da.percentile(sample,ranks).compute()
	breaks = np.convolve(centers,[0.5,0.5],mode='valid')
	centers = centers[1:-1]
	ranks = ranks[1:-1]

	return ranks, centers, breaks

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
	ranks_invlog = getInvLogRanks(n_pts,n_pts_per_bin,n_bins_per_decade)

	ranks_invlog, centers_invlog, breaks_invlog = \
		computePercentilesAndBinsFromRanks(sample,ranks_invlog)

	return ranks_invlog, centers_invlog, breaks_invlog

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
def defineLinearBins(sample,n_lin_bins=nlb,vmin=None,vmax=None):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- vmin and vmax: extremum values
	Returns:
		- centers (corresponding percentiles, or bin centers)
		- breaks (histogram bin edges)"""

	vmin,vmax = findExtremumValues(sample,vmin,vmax)
	breaks = np.linspace(vmin,vmax,n_lin_bins)
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
def compute1dDensities(sample,mode='linear',n_pts_per_bin=nppb,\
	n_bins_per_decade=nbpd,n_lin_bins=nlb,vmin=None,vmax=None):

	"""Arguments:
		- sample: 1D numpy or dask array of values
		- mode: 'linear', 'log', 'invlogQ'
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
			(used for invlogQ mode only)
		- n_bins is only used for linear mode
	Returns:
		- ranks, centers, breaks and probability densities"""

	if mode == 'linear':
		centers, breaks = defineLinearBins(sample,n_lin_bins,vmin,vmax)
		ranks = computePercentileRanksFromBins(sample,centers)
	elif mode == 'log':
		centers, breaks = defineLogBins(sample,n_bins_per_decade,vmin,vmax)
		ranks = computePercentileRanksFromBins(sample,centers)
	elif mode == 'invlogQ':
		ranks = getInvLogRanks(sample.size,n_pts_per_bin,n_bins_per_decade)
		ranks, centers, breaks = computePercentilesAndBinsFromRanks(sample,ranks)
	else:
		print("Unknown mode chosen:", mode)
		return

	if isinstance(sample,np.ndarray):
		densities, edges = np.histogram(sample,bins=breaks,density=True)
	elif isinstance(sample,da.core.Array):
		densities, edges = da.histogram(sample,bins=breaks,density=True)
	else:
		print("Unvalid data type:", type(sample))
		return

	return ranks, centers, breaks, densities

## Compute 2D bins, ranks and probability densities with transformations along axes
def compute2dDensities(sample1,sample2,mode1='linear',mode2='linear',\
	n_pts_per_bin=nppb,n_bins_per_decade=nbpd,n_lin_bins=nlb,\
	vmin1=None,vmax1=None,vmin2=None,vmax2=None):

	"""-- Yet to implement --
	Arguments: see function compute1dDensities
	Returns:
		- ranks1, centers1, breaks1, ranks2, centers2, breaks2, densities2D"""

## Extract an array of indices corresponding to a given percentile rank
def getIndicesOfPercentileRanks(values,rank,ranks,sample=None):

	"""-- Yet to implement --
	Arguments:
		- values (shaped as original data)
		- rank
		- ranks (compulsory argument in first version of the function)
		- subset (shaped as original data, if statistics is done over a subset)
	Returns:
		- 'indices' as a boolean array with same type as input 'values'"""


