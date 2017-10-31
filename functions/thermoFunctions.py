"""Module thermoFunctions

Functions to compute thermodynamic properties of the atmosphere. Adiabatic lapse
rates and temperature profiles, saturation specific humidity, 
"""


#---- Modules ----#
from math import *
import numpy as np
import numpy.ma as ma
import dask.array as da

#---- Own modules ----#
from thermoConstants import *
from daskOptions import *

#---- Functions ----#

## Air density from values of T, p and q
def airDensity(temp,pres,shum):

    """Use ideal gas law and virtual effect in linearized form
    Arguments:
        - temperature (K), pressure (Pa) and specific humidity (kg/kg) values 
        as numpy.ndarray's or dask.array's
    Returns:
        - density (kg/m3) values in the same format and type"""

    cn = getArrayType(temp)
    rho_dry = (pres)/(R_d*temp)
    virtual_coef = (eps+shum)/(eps*(1+shum))
    rho_moist = (rho_dry)/(virtual_coef)

    return rho_moist

## Saturation vapor pressure from Buck (1996)
def saturationVaporPressure(temp):

    """Argument: Temperature (K) as a numpy.ndarray or dask.array
    Returns: saturation vapor pressure (Pa) in the same format."""

    T_0 = 273.15
    
    def qvstar_numpy(temp):

        whereAreNans = np.isnan(temp)
        temp_wo_Nans = temp.copy()
        temp_wo_Nans[whereAreNans] = 0.
        # Initialize
        e_sat = np.zeros(temp.shape)
        e_sat[whereAreNans] = np.nan
        # T > 0C
        overliquid = (temp_wo_Nans > T_0)
        e_sat_overliquid = 611.21*np.exp(np.multiply(18.678-(temp-T_0)/234.5,
                                                      np.divide((temp-T_0),257.14+(temp-T_0))))
        e_sat[overliquid] = e_sat_overliquid[overliquid]
        # T < 0C 
        overice = (temp_wo_Nans < T_0)
        e_sat_overice = 611.15*np.exp(np.multiply(23.036-(temp-T_0)/333.7,
                                                   np.divide((temp-T_0),279.82+(temp-T_0))))
        e_sat[overice] = e_sat_overice[overice]

        return e_sat       # in Pa

    if temp.__class__ == np.ndarray:
        return qvstar_numpy(temp)
    elif temp.__class__ == da.core.Array:
        return temp.map_blocks(qvstar_numpy)
    elif temp.__class__ in [np.float32,float]:
        if temp > T_0:
            return 611.21*np.exp((18.678-(temp-T_0)/234.5)*(temp-T_0)/(257.14+(temp-T_0)))
        else:
            return 611.15*np.exp((23.036-(temp-T_0)/333.7)*(temp-T_0)/(279.82+(temp-T_0)))
    else:
        print("Unvalid data type:", type(temp))
        return

## Compute the saturation specific humidity based on the expressions by Buck
def saturationSpecificHumidity(temp,pres):

    """Convert from estimate of saturation vapor pressure to saturation specific
    humidity using the approximate equation qvsat ~ epsilon"""

    e_sat = saturationVaporPressure(temp)
    qvstar = (e_sat/R_v)/(pres/R_d)

    return qvstar

## Dry-adiabatic lapse rate
def dryAdiabaticLapseRate(temp,pres,spechum):

    """In pressure coordinates: Gamma_d/(rho*g) (K/Pa)."""

    cn = getArrayType(temp)
    dryGAmma_zCoord = gg/c_pd   # K/m
    rho = airDensity(temp,pres,spechum) # kg/m3
    dryGAmma_pCoord = (dryGAmma_zCoord)/(rho*gg)  # K/Pa

    return dryGAmma_pCoord

## Multiplicative factor to convert dry adiabat into moist adiabat
def moistAdiabatFactor(temp,pres):

    """Derived from the conservation of saturated moist static energy, the 
    Clausius-Clapeyron formula and the hydrostatic equation. Ignores conversion
    to ice and graupel and condensate loading."""

    cn = getArrayType(temp)
    inshape = temp.shape
    qvstar = saturationSpecificHumidity(temp,pres)
    coef_m = (1+(L_v*qvstar)/(R_d*temp))/(1+(L_v**2.*qvstar)/(c_pd*R_v*(temp**2.)))

    return coef_m

## Moist adiabatic lapse rate from condensation over liquid
def moistAdiabaticLapseRateSimple(temp,pres,spechum):

    """Returns the value of the moist adiabatic lapse rate as derived in textbooks
    from the conservation of liquid moist static energy. Convert on pressure
    coordinate by assuming hydrostaticity (K/Pa)."""

    cn = getArrayType(temp)
    Gamma_d_pCoord = dryAdiabaticLapseRate(temp,pres,spechum)  # K/Pa
    coef_m = moistAdiabatFactor(temp,pres)          # Unitless
    Gamma_m_pCoord = (Gamma_d_pCoord*coef_m) # K/Pa

    return Gamma_m_pCoord   

## Moist adiabatic temperature profile on sigma-levels from values of 
## surface temperature, and atmospheric pressure and soecific humidity
def moistAdiabatSimple(surftemp,pres,spechum,levdim=0):

    """Vertically integrate the analytic expression for the moist adiabatic 
    lapse rate from surface values (K).
    Arguments:
        - Ts (K, dimensions: [Nt,1,Nlat,Nlon])
        - p (Pa, dimensions: [Nt,Nlev,Nlat,Nlon])
        - q (kg/kg, dimensions: [Nt,Nlev,Nlat,Nlon])
    Usage with dask:
        Make sure the vertical coordinate is in dimension 0
        Make sure the chunks are the same
        Make sure there is only one chunk in the vertical
        Execute da.map_blocks(moistAdiabatSimple,Ts,p,q)
    Careful:
        When using this function with da.map_blocks, make sure the vertical
        coordinate is not subdivided into chunks; otherwise the temperature profile 
        will show zigzags.
        This function works with the vertical coordinate in first dimension. Other
        configurations haven't been tested.
    """

    cn = getArrayType(pres)

    # sh_shape = spechum.shape
    # ts_shape = surftemp.shape
    p_shape = pres.shape
    ndims = len(pres.shape)
    ind_low = [slice(None)]*ndims
    ind_low[levdim] = 0
    ind_high = ind_low.copy()
    Nlev = p_shape[levdim]

    temp = np.empty(p_shape)
    temp[ind_low] = surftemp

    for k in range(1,Nlev):
        ind_low[levdim] = k-1
        ind_high[levdim] = k
        dTdp = moistAdiabaticLapseRateSimple(temp[ind_low],
                                             pres[ind_low],
                                             spechum[ind_low])
        dp = cn.subtract(pres[ind_high],pres[ind_low])
        temp[ind_high] = cn.add(temp[ind_low],cn.multiply(dTdp,dp))

    # # Convert to dask.array # Unnecessary when usign with map_blocks
    # if pres.__class__ == da.core.Array:
    #     temp = da.from_array(temp,chunks=pres.chunks)

    return temp





## Implement Newton's method to find the zeros of a function where
