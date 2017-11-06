"""Module scalingApproximations

Functions to compute the scaling approximations for precipitation
extremes.
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import dask.array as da

#---- Own functions ----#
# currentpath = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from daskOptions import *
from thermoConstants import *
from thermoFunctions import *
from statisticalDistributions import *
    
#---- Parameters ----#

#---- Functions ----# 

## Scaling approximation from single-level omega T and p
def singleLevelScalingFromOmegaT(omega,temp,pres,efficiency=None):

    """The input arrays (np.ndarray's or dask.array's) must not have a vertical
    dimension. 'temp' and 'pres' must correspond to the same (sigma or pressure)
    level.
    Computes pointwise omega*qvstar(T,p)"""

    qvstar = saturationSpecificHumidity(temp,pres)

    if efficiency is None:
        return -omega*qvstar/gg
    else:
        return -efficiency*omega*qvstar/gg

## Scaling approximation from single-level omega and q
def singleLevelScalingFromOmegaQ(omega,spechum,efficiency=None):

    """The input arrays (np.ndarray's or dask.array's) must not have a vertical
    dimension.
    Computes pointwise omega*q."""

    if efficiency is None:
        return -omega*spechum/gg
    else:
        return -efficiency*omega*spechum/gg


## Compute vertical integral on pressure coordinates
def verticalPressureIntegral(pres,values=None,dvdp=None,levdim=0):

    """Arguments: np.ndarray's or dask.array's with identical dimensions
    Returns: @f[ \int x \frac{dp}{g}@f] from bottom to top of atmosphere.
    It is defined negative for positive values."""

    cn = getArrayType(pres)
    dp = cn.diff(pres,axis=levdim) 

    if values is None:  # If values not given, compute weight of vertical 
                        # column; normalizing factor
        return cn.nansum(dp/gg,axis=levdim)
    elif values.__class__ == da.core.Array or values.__class__ == np.ndarray:
        nlev = pres.shape[levdim]
        val_mids = np.take(values,range(nlev-1),axis=levdim) \
        - np.take(values,range(1,nlev),axis=levdim)
        val_prod = dp*val_mids
    elif values.__class__ == list:
        nlev = pres.shape[levdim]
        val_prod = dp.copy()
        for i in range(len(values)):
            val_mids = np.take(values[i],range(nlev-1),axis=levdim) \
            - np.take(values[i],range(1,nlev),axis=levdim)
            val_prod = val_prod*val_mids
    
    dvdp_prod = 1
    if dvdp is not None:
        if dvdp.__class__ == da.core.Array or dvdp.__class__ == np.ndarray:
            dvdp_prod = dvdp
        elif dvdp.__class__ == list:
            for i in range(len(dvdp)):
                dvdp_prod = dvdp_prod*dvdp[i]

    return cn.nansum(val_prod*dvdp_prod/gg,axis=levdim)

## Compute O'Gorman & Schneider's scaling




from scipy.optimize import leastsq

def computeScalingOmegaTAtAllRanks(ranks,percentiles,omega_lev,temp_lev,pres_lev,
    pr_ref,ranks_ref=None,efficiency=None):
    
    """Returns the scaling expression computed over Q-binned predictor variables.
    Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    or efficiency is not None, in which case Q_ref is not needed."""
    
    if efficiency is None:            
        efficiency = computeEfficiencyScalingOmegaT(omega_lev,temp_lev,pres_lev,
            pr_ref,ranks_ref,ranks,percentiles)
    
    pr_sc = np.array(list(map(lambda x:computeScalingOmegaTAtRank(x,ranks,percentiles,
        omega_lev,temp_lev,pres_lev,pr_ref,efficiency=efficiency),ranks)))
    return efficiency, pr_sc

def computeScalingOmegaTAtRank(rank,ranks,percentiles,omega_lev,temp_lev,pres_lev,
    pr_ref,efficiency=1):

    omega_Q = meanXAtYPercentiles(rank,ranks,percentiles,omega_lev,pr_ref)
    temp_Q = meanXAtYPercentiles(rank,ranks,percentiles,temp_lev,pr_ref)
    pres_Q = meanXAtYPercentiles(rank,ranks,percentiles,pres_lev,pr_ref)
#     pr_ref_Q = meanXAtYPercentiles(Q,Q_IL,pr_ref,pr_ref)
    
    return singleLevelScalingFromOmegaT(omega_Q,temp_Q,pres_Q,efficiency=efficiency)
        

def leastSquareCoef(pr_sc,pr_ref):
    
    guess = pr_ref[0]/pr_sc[0]
    return leastsq(lambda x:pr_ref-x*pr_sc,guess)[0][0]
    
def computeEfficiency(pr_sc,pr_ref,Q_slice=slice(None,None)):
    
    return leastSquareCoef(pr_sc[Q_slice],pr_ref[Q_slice])

def computeEfficiencyScalingOmegaT(omega_lev,temp_lev,pres_lev,pr_ref,ranks_ref,ranks,percentiles):
    
    """Compute efficiency as the tuning coefficient between percentiles of pr_ref
    and the scaling expression derived from mean values of omega, temp and pres in
    the corresponding percentile bins (Q_ref) of the pr_ref distribution."""
    
    pr_sc_zeroeff_Qs = np.array(list(map(lambda x:computeScalingOmegaTAtRank(x,ranks,
        percentiles,omega_lev,temp_lev,pres_lev,pr_ref),ranks_ref)))
    pr_ref_Qs = np.array(list(map(lambda x:meanXAtYPercentiles(x,ranks,percentiles,
        pr_ref,pr_ref),ranks_ref)))    

    return computeEfficiency(pr_sc_zeroeff_Qs,pr_ref_Qs)


def computeScalingOmegaQAtRank(rank,ranks,percentiles,omega_lev,spechum_lev,
    pr_ref,efficiency=1):

    omega_Q = meanXAtYPercentiles(rank,ranks,percentiles,omega_lev,pr_ref)
    spechum_Q = meanXAtYPercentiles(rank,ranks,percentiles,spechum_lev,pr_ref)
    
    return singleLevelScalingFromOmegaQ(omega_Q,spechum_Q,efficiency=efficiency)
  
def computeScalingOmegaQAtAllRanks(ranks,percentiles,omega_lev,spechum_lev,
    pr_ref,ranks_ref=None,efficiency=None):
    
    """Returns the scaling expression computed over Q-binned predictor variables.
    Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    or efficiency is not None, in which case Q_ref is not needed."""
    
    if efficiency is None:            
        efficiency = computeEfficiencyScalingOmegaQ(omega_lev,spechum_lev,
            pr_ref,ranks_ref,ranks,percentiles)
    
    pr_sc = np.array(list(map(lambda x:computeScalingOmegaQAtRank(x,ranks,percentiles,
        omega_lev,spechum_lev,pr_ref,efficiency=efficiency),ranks)))
    return efficiency, pr_sc

def computeEfficiencyScalingOmegaQ(omega_lev,spechum_lev,pr_ref,ranks_ref,ranks,percentiles):
    
    """Compute efficiency as the tuning coefficient between percentiles of pr_ref
    and the scaling expression derived from mean values of omega and q in
    the corresponding percentile bins (Q_ref) of the pr_ref distribution."""
    
    pr_sc_zeroeff_Qs = np.array(list(map(lambda x:computeScalingOmegaQAtRank(x,ranks,
        percentiles,omega_lev,spechum_lev,pr_ref),ranks_ref)))
    pr_ref_Qs = np.array(list(map(lambda x:meanXAtYPercentiles(x,ranks,percentiles,
        pr_ref,pr_ref),ranks_ref)))    

    return computeEfficiency(pr_sc_zeroeff_Qs,pr_ref_Qs)


