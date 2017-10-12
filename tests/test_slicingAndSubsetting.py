import sys,os
import time

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from slicingAndSubsetting import *

if __name__ == "__main__":

	print("--Test sample1d")
	x = np.random.random(size=(100,100,100))
	subset_3D = np.random.random(size=(100,100,100))
	subset_3D = np.around(subset_3D)
	subset_3D = subset_3D.astype(bool)
	subset_2D = subset_3D[0,...]
	print('- with 3D slice')
	print("# valid points =",subset_3D.sum())
	sample_1D = sample1d(x,subset_3D)
	print('sample shape:',sample_1D.shape)
	print('- with 2D slice')
	print("# valid points =",subset_2D.sum()*subset_3D.shape[0])
	sample_1D = sample1d(x,subset_2D)
	print('sample shape:',sample_1D.shape)
	print()

	print("-- Test computeTimeHorizontalMean")
	print("*With numpy arrays")
	x = np.random.random(size=(40,30,100,100))
	ind_2D = np.random.random(size=(40,100,100))
	ind_2D = np.round(ind_2D)
	start = time.time()
	x_mean = computeTimeHorizontalMean(x,ind_2D,is_3D=True)
	end = time.time()
	print("Mean of 3D data over a 2D slice",x_mean)
	print("time elapsed:", end-start)
	ind_3D = np.random.random(size=(40,30,100,100))
	ind_3D = np.round(ind_3D)
	start = time.time()
	x_mean = computeTimeHorizontalMean(x,ind_3D,is_3D=True)
	end = time.time()
	print("Mean of 3D data over a 3D slice",x_mean)
	print("time elapsed:", end-start)
	print("-")
	print("*With dask arrays")
	x_da = da.from_array(x,chunks=x.size/4)
	print(x.size)
	ind_2D_da = da.from_array(ind_2D,chunks=x.size/4)
	start = time.time()
	x_da_mean = computeTimeHorizontalMean(x_da,ind_2D_da,is_3D=True)
	end = time.time()
	print("Mean of 3D data over a 2D slice",x_da_mean.compute())
	print(x_da_mean.__class__)
	print("time elapsed:", end-start)
	print()

	sys.exit(0)