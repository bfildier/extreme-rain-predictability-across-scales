"""Module statisticalDistributions

Contains functions to compute probability densities with various transformations
on the x-axis (logarithmic, linear and inverse-logarithmic).
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import numpy.ma as ma
import dask.array as da
from scipy.interpolate import lagrange
import warnings

#---- Own functions ----#
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))
    
from daskOptions import *

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
def getInvLogRanks(n_pts,n_pts_per_bin=nppb,n_bins_per_decade=nbpd,fill_last_decade=False):

	"""Arguments:
		- n_pts: sample size
		- n_pts_per_bin: minimum number of data points per bin used to choose the maximum percentile rank
		- n_bins_per_decade: number of ranks/bins per logarithmic decade
		- fill_last_decade: True (default is False) if want to plot
		up to 99.99 or 99.999, not some weird number in the middle of a decade.
	Returns:
		- float, 1D numpy.array"""

	# k indexes bins
	n_decades = log10(n_pts/n_pts_per_bin) 		# Maximum number of decades
	dk = 1/n_bins_per_decade
	if fill_last_decade:
		k_max = ceil(n_decades) 				# Maximum bin index
	else:
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

## Get index of rank in list of ranks
def indexOfRank(rank,ranks):
    
	"""Returns the index of the closest rank in numpy.array ranks"""

	dist_to_rank = np.absolute(np.subtract(ranks,rank*np.ones(ranks.shape)))
	mindist = dist_to_rank.min()
	return np.argmax(dist_to_rank == mindist)

## Mask points that do not correspond to percentile Q of Y
def getStencilAtRank(rank,ranks,bins,Y):

	"""Returns a numpy.array of bools corresponding to percentile rank.
	It computes the corresponding rank index and bin."""

	cn = getArrayType(Y)
	i_Q = indexOfRank(rank,ranks)
	mask_Q = cn.logical_and(Y > bins[i_Q-1],Y<= bins[i_Q])

	return mask_Q

## Convert rank (float) to rank id (string)
def rankID(rank):

	return "%2.4f"%rank

## Add np.nans where Y is not defined in reference ranks axis
def adjustRanks(Y,Yranks,ranks_ref):
    
    """Assuming that all Yranks are in ranks_ref."""

    Y_adj = np.array([np.nan]*ranks_ref.size)
    for iQ in range(Yranks.size):
        iQ_ref = indexOfRank(Yranks[iQ],ranks_ref)
        Y_adj[iQ_ref] = Y[iQ]
    
    return Y_adj

## Get rank locations from rank, ranks and bins, or rank and ranks_locations
def getRankLocations(rank,Y,ranks=None,bins=None,rank_locations=None):

	rank_id = rankID(rank)
	# If had already been computed before, get it
	if rank_locations is not None:
		# print("rank_locations is not None")
		if rank_id in rank_locations.keys():
			# print("rank_id in rank_locations.keys()")
			return rank_locations[rank_id]
	
	# If had not been computed before, compute it...
	if Y.__class__ == np.ndarray:
		stencil_Q = getStencilAtRank(rank,ranks,bins,Y)
	elif Y.__class__ == da.core.Array:
		stencil_Q = da.map_blocks(lambda x: getStencilAtRank(rank,
			ranks,bins,x),Y,dtype=bool)
	# ... and store it
	rank_locations[rank_id] = stencil_Q

	return stencil_Q

## Mean of X at locations of bins of Y
def sampleSizeAtYRank(rank,Y,ranks=None,bins=None,rank_locations=None):

	# Get rank locations
	stencil_Q = getRankLocations(rank,Y,ranks,bins,rank_locations)
	# Return 0 if empty
	if stencil_Q.sum() == 0:
		return 0
	# Otherwise count
	if Y.__class__ == np.ndarray:
		return stencil_Q.sum()
	elif Y.__class__ == da.core.Array:
		return (stencil_Q.sum()).compute()

## Sample size at all ranks
def sampleSizeAtAllRanks(targetranks,Y,ranks,bins=None,rank_locations=None):

	out = np.array([np.nan]*ranks.size)
	for rank in targetranks:
		iQ = indexOfRank(rank,ranks)
		out[iQ] = sampleSizeAtYRank(rank,Y,ranks,bins,rank_locations)
	return out

## Mean of X at locations of bins of Y
def meanXAtYRank(rank,X,Y,ranks=None,bins=None,rank_locations=None):

	# Get rank locations
	stencil_Q = getRankLocations(rank,Y,ranks,bins,rank_locations)
	# Return nan if empty
	if stencil_Q.sum() == 0:
		return np.nan
	# Return nan if only nans
	if np.isnan(X[stencil_Q]).sum() == stencil_Q.sum():
		return np.nan
	# Otherwise compute nanmean
	if Y.__class__ == np.ndarray:
		return np.nanmean(X[stencil_Q])
	elif Y.__class__ == da.core.Array:
		return da.nanmean(X[stencil_Q]).compute()

## Mean of X within given percentile bins of Y
def meanXAtAllYRanks(targetranks,X,Y,ranks,bins=None,rank_locations=None):

	out = np.array([np.nan]*ranks.size)
	for rank in targetranks:
		iQ = indexOfRank(rank,ranks)
		out[iQ] = meanXAtYRank(rank,X,Y,ranks,bins,rank_locations)
	return out

## Variance of X at locations of bins of Y
def varXAtYRank(rank,X,Y,ranks=None,bins=None,rank_locations=None):

	# Get rank locations
	stencil_Q = getRankLocations(rank,Y,ranks,bins,rank_locations)
	# Return nan if empty or singleton
	if stencil_Q.sum() <= 1:
		return np.nan
	# Return nan if number of non nans is too small
	if stencil_Q.sum() - np.isnan(X[stencil_Q]).sum() <= 1:
		return np.nan
	# Otherwise compute nanvar
	if Y.__class__ == np.ndarray:
		return np.nanvar(X[stencil_Q])
	elif Y.__class__ == da.core.Array:
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			return da.nanvar(X[stencil_Q]).compute()

## Variance of X within given precentile bins of Y
def varXAtAllYRanks(targetranks,X,Y,ranks,bins=None,rank_locations=None):

	out = np.array([np.nan]*ranks.size)
	for rank in targetranks:
		iQ = indexOfRank(rank,ranks)
		out[iQ] = varXAtYRank(rank,X,Y,ranks,bins,rank_locations)
	return out

## Covariance of X1,X2 at locations of bins of Y
def covAtYRank(rank,X1,X2,Y,ranks=None,bins=None,rank_locations=None):

	cn = getArrayType(Y)
	# Define nan covariance function
	def cov(x,y):
		x_m = cn.nanmean(x)
		y_m = cn.nanmean(y)
		return cn.nanmean((x-x_m)*(y-y_m))

	# Get rank locations
	stencil_Q = getRankLocations(rank,Y,ranks,bins,rank_locations)
	# Return nan if empty or singleton
	if stencil_Q.sum() <= 1:
		return np.nan
	# Otherwise compute nancov
	if Y.__class__ == np.ndarray:
		return cov(X1[stencil_Q],X2[stencil_Q])
	elif Y.__class__ == da.core.Array:
		# with warnings.catch_warnings():
		# 	warnings.simplefilter("ignore")
		return cov(X1[stencil_Q],X2[stencil_Q]).compute()

## Covariance of X1,X2 within given percentile bins of Y
def covAtAllYRanks(targetranks,X1,X2,Y,ranks,bins=None,rank_locations=None):

	out = np.array([np.nan]*ranks.size)
	for rank in targetranks:
		iQ = indexOfRank(rank,ranks)
		out[iQ] = covAtYRank(rank,X1,X2,Y,ranks,bins,rank_locations)
	return out

## Percentiles of X within percentile bins of Y
def XPercentilesAtYRank(rank_X,X,ranks_Y,Y,ranks_X=None,bins_X=None,
	rank_locations_X=None):

	cn = getArrayType(X)

	# Get rank locations
	stencil_Q = getRankLocations(rank_X,Y,ranks_X,bins_X,rank_locations_X)
	# Compute percentile
	X_at_rank = X[stencil_Q] if cn == np else X[stencil_Q].compute()
	if X_at_rank.size == 0 :
		out = np.array([np.nan]*len(ranks_Y))
	else:
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			out = np.percentile(X_at_rank,ranks_Y)
	return out

## Percentiles of X within given percentile bins of Y
def XPercentilesAtAllYRanks(targetranks_X,X,ranks_Y,Y,ranks_X,bins_X=None,
	rank_locations_X=None):

	out = np.array([[np.nan]*len(ranks_Y)]*ranks_X.size)
	for rank_X in targetranks_X:
		iQ = indexOfRank(rank_X,ranks_X)
		out[iQ] = XPercentilesAtYRank(rank_X,X,ranks_Y,Y,ranks_X,bins_X,
			rank_locations_X)
	return out

