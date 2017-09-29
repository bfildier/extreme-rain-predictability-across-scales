"""@package docstring
Module statisticalDistributions

Contains functions to compute probability densities with various transformations
on the x-axis (logarithmic, linear and inverse-logarithmic).
"""

###--- Modules ---###
import sys,os
from math import *
import numpy as np
import dask.array as da

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from defineDaskParameters import *
    
def getInvLogRanks(n_pts,n_pts_per_bin=4,n_bins_per_decade=10):

	"""Given a sample size n_pts and a minimum number of data points per bin,
	returns a 1D-array or percentile ranks regularly spaced on an inverse-
	logarithmic axis"""

	n_max = int(log10(n_pts/n_pts_per_bin))
	di = 1/n_bins_per_decade
	max_rank = (1.-10**(-n_max))*100
	scale_invlog = np.arange(0,n_max+di,di)
	ranks_invlog = np.subtract(np.ones(scale_invlog.size),np.power(10,-scale_invlog))*100

	return max_rank, ranks_invlog

def getSample1D(values,subset=None):

	"""Return a 1D-array containing the relavant data points.
	It assumes that values has dimensions (time,lat,lon) 
	and subset either (lat,lon) or (time,lat,lon)."""

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

def computePercentilesAndBinsFromRanks(sample,ranks):

	"""From a sequence of percentile ranks and the data sample, compute and
	returns the corresponding percentiles of the distributions and the bin edges."""

	centers = np.percentile(sample,ranks)
	breaks = (centers[1:]+centers[:-1])/2
	centers = centers[1:-1]
	ranks = ranks[1:-1]

	return ranks, centers, breaks

def computePercentilesAndBinsFromRanks_da(sample,ranks):

	"""From a sequence of percentile ranks and the data sample, compute and
	returns the corresponding percentiles of the distributions and the bin edges."""

	centers = da.percentile(sample,ranks).compute()
	breaks = (centers[1:]+centers[:-1])/2
	centers = centers[1:-1]
	ranks = ranks[1:-1]

	return ranks, centers, breaks

## Define inverse-logarithmic percentiles and bin edges
## Return: maximum percentile, percentiles (independent from values)
##     and centers_invlog (information about values)
## n_pts_per_bin is the minimum number of values per bin, used to compute highest percentile
def defineInvLogPercentiles(sample,subset=None,n_pts_per_bin=4,n_bins_per_decade=10):

	"""Compute percentiles of the distribution on a inverse-logarithmic axis"""
	
	n_pts = sample.size
	max_rank, ranks_invlog = getInvLogRanks(n_pts,n_pts_per_bin,n_bins_per_decade)

	ranks_invlog, centers_invlog, breaks_invlog = \
		computePercentilesAndBinsFromRanks(sample,ranks_invlog)

	return max_rank, ranks_invlog, centers_invlog, breaks_invlog


	
