## Modules
from math import *
import numpy as np
import numpy.ma as ma

## Own modules
# from physicalConstants import *
from thermoConstants import *
g = gg

## Air density profile based on profiles for T, p and q
## Use ideal gas law and virtual effect in linearized form
def airDensity(temp,pres,shum):
    rho_dry = np.divide(pres,R_d*temp)
    return np.divide(rho_dry,(1+eps*shum))

# Compute the saturation vapor pressure from expressions by Buck (1996)
def saturationVaporPressure(temp):
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

# Compute the saturation specific humidity based on the expressions by Buck
def saturationSpecificHumidity(temp,pres):
    e_sat = saturationVaporPressure(temp)
    return np.divide(e_sat/R_v,pres/R_d)