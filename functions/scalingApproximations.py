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
        return omega*qvstar
    else:
        return efficiency*omega*qvstar

## Scaling approximation from single-level omega and q
def singleLevelScalingFromOmegaQ(spechum,omega,efficiency=None):

    """The input arrays (np.ndarray's or dask.array's) must not have a vertical
    dimension.
    Computes pointwise omega*q."""

    if efficiency is None:
        return omega*spechum
    else:
        return efficiency*omega*spechum


## Compute vertical integral on pressure coordinates
def verticalPressureIntegral(pres,values=None,levdim=0):

    """Arguments: np.ndarray's or dask.array's with identical dimensions
    Returns: @f[ \int x \frac{dp}{g}@f] from bottom to top of atmosphere.
    It is defined negative for positive values."""

    cn = getArrayType(pres)
    dp = cn.diff(pres,axis=levdim) 

    if values is None:  # If values not given, compute weight of vertical 
                        # column; normalizing factor
        return cn.nansum(dp/gg,axis=levdim)
    else:   # Compute integral
        nlev = pres.shape[levdim]
        val_mids = np.take(values,range(nlev-1),axis=levdim) \
        - np.take(values,range(1,nlev),axis=levdim)
        return cn.nansum(val_mids*dp/gg,axis=levdim)

