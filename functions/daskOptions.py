"""Module daskOptions

Define parameters related to the use of dask arrays that are common across modules
"""

#---- Modules ----#
import numpy as np
import dask.array as da

## Default compute action for dask arrays
da_compute_default = False
## Size of dask array chunks
chunks = 100000

#---- Functions ----#

## Get module corresponding to array type
def getArrayType(array):

	if array.__class__.__bases__[0] is np.ndarray or array.__class__ is np.ndarray:
		cn = np # stands for 'class name'
	elif array.__class__ == da.core.Array:
		cn = da
	else:
		print("[Error in daskOptions.getArrayType] Unvalid data type:", type(array))
		return 

	return cn