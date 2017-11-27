"""Module scalingApproximations

Functions to compute the scaling approximations for precipitation
extremes.
"""

#---- Modules -----#
import sys,os
from math import *
import numpy as np
import dask.array as da
from scipy.optimize import leastsq,curve_fit

#---- Own functions ----#
# currentpath = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from daskOptions import *
from thermoConstants import *
from thermoFunctions import *
from statisticalDistributions import *
    
#---- Parameters ----#

#---- Functions ----# 

##-- Single-level scaling calculations --#

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


def computeScalingOmegaTAtAllRanks(ranks,omega_lev,temp_lev,pres_lev,pr_ref,
    ranks_ref=None,percentiles_ref=None,efficiency=None,bins=None,
    rank_locations=None):
    
    """Returns the scaling expression computed over Q-binned predictor variables.
    Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    or efficiency is not None, in which case Q_ref is not needed."""
    
    if efficiency is None:            
        efficiency = computeEfficiencyScalingOmegaT(omega_lev,temp_lev,pres_lev,
            pr_ref,ranks_ref,percentiles_ref,ranks,bins,rank_locations)
    
    pr_sc = np.array(list(map(lambda x:computeScalingOmegaTAtRank(x,omega_lev,
        temp_lev,pres_lev,pr_ref,efficiency,ranks,bins,rank_locations),
        ranks)))

    return efficiency, pr_sc

def computeScalingOmegaTAtRank(rank,omega_lev,temp_lev,pres_lev,pr_ref,
    efficiency=1,ranks=None,bins=None,rank_locations=None):

    omega_Q = meanXAtYRank(rank,omega_lev,pr_ref,ranks,bins,rank_locations)
    temp_Q = meanXAtYRank(rank,temp_lev,pr_ref,ranks,bins,rank_locations)
    pres_Q = meanXAtYRank(rank,pres_lev,pr_ref,ranks,bins,rank_locations)
    
    return singleLevelScalingFromOmegaT(omega_Q,temp_Q,pres_Q,efficiency=efficiency)
        

def leastSquareCoef(pr_sc,pr_ref):
    
    guess = pr_ref[0]/pr_sc[0]
    def fun(p):
        diff_p = pr_ref-p*pr_sc
        diff_p[np.isnan(diff_p)] = 0
        return diff_p

    return leastsq(fun,guess)[0][0]
    
def computeEfficiency(pr_sc,pr_ref,Q_slice=slice(None,None)):
    
    return leastSquareCoef(pr_sc[Q_slice],pr_ref[Q_slice])

def computeEfficiencyScalingOmegaT(omega_lev,temp_lev,pres_lev,pr_ref,ranks_ref,
    percentiles_ref,ranks=None,bins=None,rank_locations=None):
    
    """Compute efficiency as the tuning coefficient between bins of pr_ref
    and the scaling expression derived from mean values of omega, temp and pres in
    the corresponding percentile bins (Q_ref) of the pr_ref distribution."""
    
    pr_sc_zeroeff_Qs = np.array(list(map(lambda x:computeScalingOmegaTAtRank(x,
        omega_lev,temp_lev,pres_lev,pr_ref,1,ranks,bins,rank_locations),
        ranks_ref)))
    # pr_ref_Qs = np.array(list(map(lambda x:meanXAtYRank(x,pr_ref,pr_ref,
        # ranks,bins,rank_locations),ranks_ref)))    
    pr_ref_Qs = percentiles_ref

    return computeEfficiency(pr_sc_zeroeff_Qs,pr_ref_Qs)

def computeScalingOmegaQAtRank(rank,omega_lev,spechum_lev,pr_ref,efficiency=1,
    ranks=None,bins=None,rank_locations=None):

    omega_Q = meanXAtYRank(rank,omega_lev,pr_ref,ranks,bins,rank_locations)
    spechum_Q = meanXAtYRank(rank,spechum_lev,pr_ref,ranks,bins,rank_locations)
    
    return singleLevelScalingFromOmegaQ(omega_Q,spechum_Q,efficiency=efficiency)
  
def computeScalingOmegaQAtAllRanks(ranks,omega_lev,spechum_lev,pr_ref,
    ranks_ref=None,percentiles_ref=None,efficiency=None,bins=None,
    rank_locations=None):
    
    """Returns the scaling expression computed over Q-binned predictor variables.
    Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    or efficiency is not None, in which case Q_ref is not needed."""
    
    if efficiency is None:            
        efficiency = computeEfficiencyScalingOmegaQ(omega_lev,spechum_lev,
            pr_ref,ranks_ref,percentiles_ref,ranks,bins,rank_locations)
    
    pr_sc = np.array(list(map(lambda x:computeScalingOmegaQAtRank(x,omega_lev,
        spechum_lev,pr_ref,efficiency,ranks,bins,rank_locations),ranks)))

    return efficiency, pr_sc

def computeEfficiencyScalingOmegaQ(omega_lev,spechum_lev,pr_ref,ranks_ref,
    percentiles_ref,ranks=None,bins=None,rank_locations=None):
    
    """Compute efficiency as the tuning coefficient between bins of pr_ref
    and the scaling expression derived from mean values of omega and q in
    the corresponding percentile bins (Q_ref) of the pr_ref distribution."""
    
    pr_sc_zeroeff_Qs = np.array(list(map(lambda x:computeScalingOmegaQAtRank(x,
        omega_lev,spechum_lev,pr_ref,1,ranks,bins,rank_locations),ranks_ref)))
    # pr_ref_Qs = np.array(list(map(lambda x:meanXAtYRank(x,pr_ref,pr_ref,
        # ranks,bins,rank_locations),ranks_ref)))    
    pr_ref_Qs = percentiles_ref

    return computeEfficiency(pr_sc_zeroeff_Qs,pr_ref_Qs)


##-- O'Gorman&Schneider scaling --##

## Vertical integral on pressure coordinates
def verticalPressureIntegral(pres,values=None,dvdp=None,levdim=0):

    """Arguments: np.ndarray's or dask.array's with identical dimensions
    Returns: @f[ \int x \frac{dp}{g}@f] from bottom to top of atmosphere.
    It is defined negative for positive values.
    """

    cn = getArrayType(pres)
    dp = cn.diff(pres,axis=levdim) 

    val_prod = dp.copy()
    if values is not None:
        if values.__class__ == da.core.Array or values.__class__ == np.ndarray:
            nlev = pres.shape[levdim]
            val_mids = (np.take(values,range(nlev-1),axis=levdim) \
            + np.take(values,range(1,nlev),axis=levdim))/2
            val_prod = dp*val_mids
        if values.__class__ == list:
            nlev = pres.shape[levdim]
            val_prod = dp.copy()
            for i in range(len(values)):
                val_mids = (np.take(values[i],range(nlev-1),axis=levdim) \
                + np.take(values[i],range(1,nlev),axis=levdim))/2
                val_prod = val_prod*val_mids
    
    dvdp_prod = 1
    if dvdp is not None:
        if dvdp.__class__ == da.core.Array or dvdp.__class__ == np.ndarray:
            dvdp_prod = dvdp
        elif dvdp.__class__ == list:
            for i in range(len(dvdp)):
                dvdp_prod = dvdp_prod*dvdp[i]

    return cn.nansum(val_prod*dvdp_prod/gg,axis=levdim)

## Find cold-point tropopause level
def tropopauseIndex(temp,levdim=0):
    
    cn = getArrayType(temp)

    return cn.argmin(temp,axis=levdim)
    # (or return when temperature first stops decreasing)
    # return temp1D.size - cn.argmax(cn.diff(cn.flipud(temp1D)) > 0)

## Find bottom index
def bottomIndex(pres,levdim=0):

    cn = getArrayType(pres)

    return cn.argmax(pres,axis=levdim)

## Crop profiles between surface and tropopause
def cropProfiles(pres,temp,values=None,levdim=0):

    """Only developed for numpy arrays"""

    cn = getArrayType(pres)

    # Replace 1D-data values v with nans outside slice s
    def cropVector(v,s):
        v_c = v.copy()
        v_c[:] = np.nan
        v_c[s] = v[s]
        return v_c
    # Apply cropVector function to all dimensions
    def cropArray(v,i_bot,i_trop,levdim):
        v_c = np.moveaxis(v,levdim,-1).copy()
        vshape = v_c.shape[:-1]
        for ind in cn.ndindex(vshape):
            i_min = min(i_bot[ind],i_trop[ind])
            i_max = max(i_bot[ind],i_trop[ind])
            s = slice(i_min,i_max+1)
            v_c[ind] = cropVector(v_c[ind],s)
        return np.moveaxis(v_c,-1,levdim)

    i_bot = bottomIndex(pres,levdim=levdim)
    i_trop = tropopauseIndex(temp,levdim=levdim)
    # Crop pressure and temperature arrays
    pres_c = cropArray(pres,i_bot,i_trop,levdim)
    temp_c = cropArray(temp,i_bot,i_trop,levdim)
    # Crop all other values arrays
    if values.__class__ == np.ndarray or values.__class__ == da.core.Array:
        values_c = cropArray(values,i_bot,i_trop,levdim)
        return pres_c, temp_c, values_c
    elif values.__class__ == list:
        values_c = []
        for vals in values:
            values_c.append(cropArray(vals,i_bot,i_trop,levdim))
        return (pres_c,temp_c)+tuple(values_c)

## Compute O'Gorman & Schneider's scaling
def scalingOGS09(omega,temp,pres,efficiency=1,temp_type='profile',parameter=1,
    levdim=0):
    
    """Has not been tested for dask arrays."""

    cn = getArrayType(pres)

    if temp_type == 'profile':
        temp_profile = temp
    else:
        print("wrong type of temperature profile requested. Code this option.")
    # elif temp_type == 'parametric':
    #     p_ref = cn.take(pres,-1,axis=levdim)
    #     temp_profile = 

    pres_c, temp_c, omega_c = cropProfiles(pres,temp_profile,omega,levdim=levdim)
    qvstar_c = saturationSpecificHumidity(temp_c,pres_c)
    dqvstar_dp_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(pres_c,axis=levdim)

    return -verticalPressureIntegral(pres_c,
                                     values=omega_c,
                                     dvdp=dqvstar_dp_c,
                                     levdim=levdim)

## Compute OGS09 scaling within a given percentile bin
def computeScalingOGS09AtRank(rank,omega,temp,pres,pr_ref,temp_type='profile',
    parameter=1,efficiency=1,ranks=None,bins=None,rank_locations=None):
    
    omega_Q = meanXProfileAtYRank(rank,omega,pr_ref,ranks,bins,rank_locations)
    temp_Q = meanXProfileAtYRank(rank,temp,pr_ref,ranks,bins,rank_locations)
    pres_Q = meanXProfileAtYRank(rank,pres,pr_ref,ranks,bins,rank_locations)
    
    if pres_Q is np.nan or temp_Q is np.nan or omega_Q is np.nan:
        return np.nan

    return scalingOGS09(omega_Q,temp_Q,pres_Q,efficiency,temp_type,parameter,
                        levdim=0)

## Efficiency of OGS09 scaling
def computeEfficiencyScalingOGS09(omega,temp,pres,pr_ref,ranks_ref,
    percentiles_ref=None,temp_type='profile',parameter=1,ranks=None,bins=None,
    rank_locations=None):
    
    """Compute efficiency as the tuning coefficient between bins of pr_ref
    and the scaling expression derived from mean profiles of omega, T and p in
    the corresponding percentile bins (ranks_ref) of the pr_ref distribution."""
    
    pr_sc_zeroeff_Qs = np.array(list(map(lambda x:computeScalingOGS09AtRank(x,
        omega,temp,pres,pr_ref,temp_type,parameter,1,ranks,bins,
        rank_locations),ranks_ref)))
    # pr_ref_Qs = np.array(list(map(lambda x:meanXAtYRank(x,pr_ref,pr_ref,
        # ranks,bins,rank_locations),ranks_ref)))
    pr_ref_Qs = percentiles_ref

    return computeEfficiency(pr_sc_zeroeff_Qs,pr_ref_Qs)

## Scaling at all ranks
def computeScalingOGS09AtAllRanks(ranks,omega,temp,pres,pr_ref,
    temp_type='profile',parameter=1,ranks_ref=None,percentiles_ref=None,
    efficiency=None,bins=None,rank_locations=None):
    
    """Returns the scaling expression computed over Q-binned predictor variables.
    Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    or efficiency is not None, in which case Q_ref is not needed."""
    
    if efficiency is None:            
        efficiency = computeEfficiencyScalingOGS09(omega,temp,pres,pr_ref,
            ranks_ref,percentiles_ref,temp_type,parameter,ranks,bins,rank_locations)
    
    pr_sc = np.array(list(map(lambda x:computeScalingOGS09AtRank(x,omega,temp,
        pres,pr_ref,temp_type,parameter,efficiency,ranks,bins,rank_locations),ranks)))

    return efficiency, pr_sc


#---- New scaling approximation ----#

## Extension to O'Gorman & Schneider scaling approximation
def scalingRH(omega,temp,pres,relhum,fracarea_boost=1,entrainment=1,
    levdim=1,temp_type='profile',parameter=1):
    
    """Adding a term that depends on relative humidity, and interpret the first
    coefficient as the ratio between effective updraft velocity and GCM-scale 
    vertical velocity."""
    
    cn = getArrayType(omega)

    if temp_type == 'profile':
        temp = temp

    qvstar = saturationSpecificHumidity(temp,pres)
    spechum = relhum*qvstar
    rho = airDensity(temp,pres,spechum)
    
    pres_c, temp_c, omega_c, qvstar_c, rho_c, relhum_c = cropProfiles(pres,temp,
        [omega,qvstar,rho,relhum],levdim=levdim)
    dqvstar_dp_c = cn.diff(qvstar_c,axis=levdim)/cn.diff(pres_c,axis=levdim)
    
    pr_sc = fracarea_boost * (- verticalPressureIntegral(pres_c,omega_c,
        dqvstar_dp_c,levdim=levdim) +\
    entrainment*verticalPressureIntegral(pres_c,[omega_c,qvstar_c/rho_c/gg,
        (1-relhum_c/100)],levdim=levdim))
#     pr_sc = 1/(2*fracarea_eff-1)*(- verticalPressureIntegral(pres_c,values=omega_c,dvdp=dqvstar_dp_c,levdim=levdim) +\
#     entrainment*verticalPressureIntegral(pres_c,values=[omega_c,qvstar_c/rho_c/gg,(1-relhum_c/100)],levdim=levdim))
    
    return pr_sc

## Compute parameters of the extended OGS09 scaling
def computeParametersScalingRH(omega,temp,pres,relhum,pr,ranks_ref,
    temp_type='profile',parameter=1,ranks=None,bins=None,
    rank_locations=None):
    
    cn = getArrayType(omega)

    #-- Create optimization function --##
    
    def func(x,a,b):
        s = []
        for i in range(4):
            s.append(slice(i*nlev,(i+1)*nlev))
        return scalingRH(x[s[0]],x[s[1]],x[s[2]],x[s[3]],fracarea_boost=a,
                         entrainment=b,levdim=0,temp_type=temp_type,parameter=parameter)
    
    #-- Get variables at reference locations --#
    
    varnames_for_scRH = ('omega','temp','pres','relhum','pr')

    # Initialize list for each variable
    for varname in varnames_for_scRH:
        # setattr(thismodule,"%s_ref_list"%varname,[])
        locals()["%s_ref_list"%varname] = []
    # fill list with sample variables at each percentile
    for rank in ranks_ref:
        stencil_Q = getRankLocations(rank,pr,ranks,bins,rank_locations)
        for varname in varnames_for_scRH:
            var = locals()[varname]
            # var_list = getattr(thismodule,"%s_ref_list"%varname)
            var_list = locals()["%s_ref_list"%varname]
            var_list.append(sampleFlattened(var,stencil_Q))
    # concat arrays from list
    for varname in varnames_for_scRH:
        # var_list = getattr(thismodule,"%s_ref_list"%varname)
        var_list = locals()["%s_ref_list"%varname]
        # setattr(thismodule,"%s_ref"%varname,cn.hstack(var_list))
        locals()["%s_ref"%varname] = cn.hstack(var_list)
    
    #-- Fit curve to precipitation values --#

    nlev = locals()['omega_ref'].shape[0]
    ndata = locals()['omega_ref'].shape[1]
    xdata = np.vstack([locals()['omega_ref'],
                       locals()['temp_ref'],
                       locals()['pres_ref'],
                       locals()['relhum_ref']])
    ydata = locals()['pr_ref']
    p,c = curve_fit(func,xdata,ydata)
    
    return p,c

## Compute OGS09 scaling within a given percentile bin
def computeScalingRHAtRank(rank,omega,temp,pres,relhum,pr,temp_type='profile',
    parameter=1,fracarea_boost=1,entrainment=1,
    ranks=None,bins=None,rank_locations=None):

    omega_Q = meanXProfileAtYRank(rank,omega,pr,ranks,bins,rank_locations)
    temp_Q = meanXProfileAtYRank(rank,temp,pr,ranks,bins,rank_locations)
    pres_Q = meanXProfileAtYRank(rank,pres,pr,ranks,bins,rank_locations)
    relhum_Q = meanXProfileAtYRank(rank,relhum,pr,ranks,bins,rank_locations)

    if pres_Q is np.nan or temp_Q is np.nan or omega_Q is np.nan:
        return np.nan

    return scalingRH(omega_Q,temp_Q,pres_Q,relhum_Q,fracarea_boost=fracarea_boost,
        entrainment=entrainment,levdim=0,temp_type=temp_type,parameter=parameter)

## Scaling at all ranks
def computeScalingRHAtAllRanks(ranks,omega,temp,pres,relhum,pr,temp_type='profile',
    parameter=1,fracarea_boost=None,entrainment=None,
    ranks_ref=None,bins=None,rank_locations=None):
    
    # """Returns the scaling expression computed over Q-binned predictor variables.
    # Either efficiency is None, in which case need to compute efficiency for Q_ref range,
    # or efficiency is not None, in which case Q_ref is not needed."""
    
    if fracarea_boost is None or entrainment is None:
        p,c = computeParametersScalingRH(omega,temp,pres,relhum,pr,ranks_ref,
            temp_type=temp_type,parameter=parameter,ranks=ranks,bins=bins,
            rank_locations=rank_locations)
        fracarea_boost = p[0]
        entrainment = p[1]
    
    pr_sc = np.array(list(map(lambda x:computeScalingRHAtRank(x,omega,temp,pres,
        relhum,pr,temp_type=temp_type,parameter=parameter,
        fracarea_boost=fracarea_boost,entrainment=entrainment,ranks=ranks,
        bins=bins,rank_locations=rank_locations),ranks)))

    return fracarea_boost, entrainment, pr_sc

