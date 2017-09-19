import sys,os
import numpy as np

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

from environmentAndDirectories import *
from importingData import *

if __name__ == "__main__":

	print("--Test getInputDirectories")
	dataroot = '/Users/bfildier/Data'
	compset = 'FSPCAMm_AMIP'
	experiment = 'piControl'
	dates = ('1850050100000','1850050200000')
	print("dataroot   :",dataroot)
	print("compset    :",compset)
	print("experiment :",experiment)
	print("dates      :",dates)

	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results,\
	 inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	print(inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results, \
		inputdir_fx)
	print()

	print("--Test getProcessedFiles")
	inputfiles = getProcessedFiles('PRECT',inputdir_processed_1hr)
	print("Hourly data:", inputfiles)
	print()

	print("--Test getProcessedValues")
	print("* With required arguments only:")
	values = getProcessedValues('PRECT', inputdir_processed_day)
	print("values.shape =", values.shape)
	print("* With specified inputfiles (bypass argument 'inputdir'):")
	values = getProcessedValues('PRECT', inputdir_processed_day, inputfiles = inputfiles)
	print("values.shape =", values.shape)
	print()

	print("--Test getSimulationValues")
	values = getSimulationValues('PRECT',inputdir,inputfiles=None,dates=dates,handle='h0')
	if values is not None:
		print("values.shape =", values.shape)
	print()

	print("--Test getValues")
	for varid in 'PRECT','I50':
		print("* With %s:"%varid)
		values = getValues(varid,dataroot,compset,subset='tropics',\
			experiment=experiment,time_stride='1hr',time_type='A',dates=dates)
		if values is not None:
			print("values.shape =",values.shape)
	print()

	print("--Test getPressureCoordinateFunction")
	lev_file = os.path.join(inputdir_fx,"lev_fx_CESM111-SPCAM20_allExperiments_r0i0p0.nc")
	p0_scalar = 100000
	p0_array = p0_scalar*np.ones((2,3))
	computeP = getPressureCoordinateFunction(lev_file)
	print("* Scalar:")
	print("1D pressure values:\n",computeP(p0_scalar))
	print("* numpy array:")
	print("pressure.shape =",computeP(p0_array).shape)
	print()

	sys.exit(0)