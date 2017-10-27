import sys,os
import numpy as np
import dask.array as da
import matplotlib.pyplot as plt
import time

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from scalingApproximations import *

if __name__ == "__main__":

	print("-- Test verticalPressureIntegral --")
	print()
	print("- 1D - ")
	pres = np.linspace(100000,0,30)
	print(pres)
	pres_da = da.from_array(pres,chunks=15)
	print("Mass of atmospheric column is",verticalPressureIntegral(pres_da).compute(),"kg/m2")
	print("Pressure-weighted pressure:", verticalPressureIntegral(pres_da,pres_da).compute(),"kg2/m3/s2")
	print()
	print("- 2D, lev last -")
	x = np.linspace(0,100,3)
	# x_da = da.from_array(x,chunks=5)
	pres_2D, x_2D = np.meshgrid(pres,x)
	print("shape:",pres_2D.shape)
	levdim = 1
	pres_2D_da = da.from_array(pres_2D,chunks=(5,15))
	x_2D_da = da.from_array(x_2D,chunks=(5,15))
	print("Mass of atmospheric column is",verticalPressureIntegral(pres_2D_da,
		levdim=levdim).compute(),"kg/m2")
	print("Pressure-weighted pressure:", verticalPressureIntegral(pres_2D_da,
		pres_2D_da,levdim=levdim).compute(),"kg2/m3/s2")
	print("Pressure-weighted pressure^2, args as list:",
		verticalPressureIntegral(pres_2D_da,
		[pres_2D_da,pres_2D_da],levdim=levdim).compute(),"kg3/m4/s4")
	print()
	print("- 2D, lev first -")
	x = np.linspace(0,100,3)
	# x_da = da.from_array(x,chunks=5)
	x_2D, pres_2D = np.meshgrid(x,pres)
	print("shape:",pres_2D.shape)
	levdim = 0
	pres_2D_da = da.from_array(pres_2D,chunks=(5,15))
	x_2D_da = da.from_array(x_2D,chunks=(5,15))
	print("Mass of atmospheric column is",verticalPressureIntegral(pres_2D_da,
		levdim=levdim).compute(),"kg/m2")
	print("Pressure-weighted pressure:", verticalPressureIntegral(pres_2D_da,
		pres_2D_da,levdim=levdim).compute(),"kg2/m3/s2")
	print("Pressure-weighted pressure^2, args as list:",
		verticalPressureIntegral(pres_2D_da,
		[pres_2D_da,pres_2D_da],levdim=levdim).compute(),"kg3/m4/s4")
	print()

	
	sys.exit(0)