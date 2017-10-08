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

## Air density from values of T, p and q
def airDensity(temp,pres,shum):

    """Use ideal gas law and virtual effect in linearized form
    Arguments:
        - temperature (K), pressure (Pa) and specific humidity (kg/kg) values 
        as numpy.ndarray's or dask.array's
    Returns:
        - density (kg/m3) values in the same format and type"""

    print("ADAPT THIS FUNCTION FOR DASK")

    rho_dry = np.divide(pres,R_d*temp)
    virtual_coef = np.divide(eps+shum,eps*(1+shum))
    return np.divide(rho_dry,virtual_coef)

## Saturation vapor pressure from Buck (1996)
def saturationVaporPressure(temp):

    """Argument: Temperature (K) as a numpy.ndarray or dask.array
    Returns: saturation vapor pressure (Pa) in the same format."""

    print("ADAPT THIS FUNCTION FOR DASK")

    T_0 = 273.15
    whereAreNans = np.isnan(temp)
    temp_wo_Nans = temp.copy()
    temp_wo_Nans[whereAreNans] = 0.
    # Initialize
    e_sat = np.zeros(temp.shape)
    e_sat[whereAreNans] = np.nan
    # T > 0C
    overliquid = np.array((temp_wo_Nans > T_0),dtype=bool)
    e_sat_overliquid = 0.61121*np.exp(np.multiply(18.678-(temp-T_0)/234.5,
                                                  np.divide((temp-T_0),257.14+(temp-T_0))))
    e_sat[overliquid] = e_sat_overliquid[overliquid]
    # T < 0C 
    overice = np.array((temp_wo_Nans < T_0),dtype=bool)
    e_sat_overice = 0.61115*np.exp(np.multiply(23.036-(temp-T_0)/333.7,
                                               np.divide((temp-T_0),279.82+(temp-T_0))))
    e_sat[overice] = e_sat_overice[overice]
    
    return e_sat*1000       # in Pa

## Compute the saturation specific humidity based on the expressions by Buck
def saturationSpecificHumidity(temp,pres):

    """Convert from estimate of saturation vapor pressure to saturation specific
    humidity using the approximate equation qvsat ~ epsilon"""

    print("ADAPT THIS FUNCTION FOR DASK")

    e_sat = saturationVaporPressure(temp)
    return np.divide(e_sat/R_v,pres/R_d)

## Dry-adiabatic lapse rate
def dryAdiabaticLapseRate(temp,pres,spechum):

    """In pressure coordinates: Gamma_d/(rho*g) (K/Pa)."""

    dryGAmma_zCoord = gg/c_pd   # K/m
    rho = airDensity(temp,pres,spechum) # kg/m3
    dryGAmma_pCoord = dryGAmma_zCoord/(rho*gg)  # K/Pa

    return dryGAmma_pCoord

## Multiplicative factor to convert dry adiabat into moist adiabat
def moistAdiabatFactor(temp,pres):

    """Derived from the conservation of saturated moist static energy, the 
    Clausius-Clapeyron formula and the hydrostatic equation. Ignores conversion
    to ice and graupel and condensate loading."""

    print("ADAPT THIS FUNCTION FOR DASK")

    inshape = temp.shape
    qvstar = saturationSpecificHumidity(temp,pres)

    coef_m = np.divide(np.ones(inshape)+np.divide(L_v*qvstar,
                                                 R_d*temp),
                      np.ones(inshape)+np.divide(L_v**2.*qvstar,
                                                 c_pd*R_v*np.power(temp,2.)))

    return coef_m

## Moist adiabatic lapse rate from condensation over liquid
def moistAdiabaticLapseRateLiquid(temp,pres,spechum):

    """Returns the value of the moist adiabatic lapse rate as derived in textbooks
    from the conservation of liquid moist static energy. Convert on pressure
    coordinate by assuming hydrostaticity (K/Pa)."""

    print("ADAPT THIS FUNCTION FOR DASK")

    Gamma_d_pCoord = dryAdiabaticLapseRate(temp,pres,spechum)  # K/Pa
    coef_m = moistAdiabatFactor(temp,pres)          # Unitless

    return Gamma_d_pCoord*coef_m    # K/Pa

## Moist adiabatic temperature profile on sigma-levels from values of 
## surface temperature, and atmospheric pressure and soecific humidity
def moistAdiabat(surftemp,pres,spechum):

    """Vertically integrate the analytic expression for the moist adiabatic 
    lapse rate from surface values (K)."""


## Implement Newton's method to find the zeros of a function where