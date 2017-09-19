###--- Modules ---###

import glob
import sys,os
import string
import numpy as np
from datetime import date, time, datetime, timedelta
from netCDF4 import Dataset

## Own functions
# currentpath = os.path.dirname(os.path.realpath(__file__))
# 		

from environmentAndDirectories import *

###--- Functions ---###

## inputfiles is a list of inputfiles to open simultaneously.
## dates is a pair of date strings in the form 'YYYYMMDD'
## - if both arguments 'inputfiles' and 'dates' are NA, 
##   use all files $inputdir/$varid_*
## - if inputfiles is given, use this list in priority
## - if inputfiles is not given but dates is given, use corresponding files
## between those dates
#
# 1. Get inputfiles in correct format
def getInputfiles(varid,inputdir,inputfiles=None,dates=None):

	"""Based on a variable ID and input directory, output a list of corresponding 
	filenames for preprocessed files.
	 - filenames have to start with $varid
	 - filenames must contain the time as .YYYYMMDD.nc"""

	## Get the list of input files
	if inputfiles == None:
		if dates == None:    # Get all matching files
			inputfiles = glob.glob(os.path.join(inputdir,varid.replace('_','-')+"_*.nc"))
		else:    # Filter between dates
			dt_bnds = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8])) for d in dates]
			inputfiles = []
			for file in glob.glob(os.path.join(inputdir,varid.replace('_','-')+"_*.nc")):
				filename = os.path.basename(file)
				dt_info = filename.split('.')[0].split('_')[-1].split('-')
				dt = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8])) for d in dt_info]
				if dt[0] >= dt_bnds[0] and dt[1] <= dt_bnds[1]:
					inputfiles.append(file)
	else:
		## Append dirname to all files if necessary
		inputfiles = [os.path.join(inputdir,f) if inputdir not in f else f for f in inputfiles]
	inputfiles.sort()
	if len(inputfiles) == 0:
		print(inputdir, varid, dates)
	return inputfiles
#
# 2. Get values for processed data ($dataroot/preprocessed/$case/$freq/*)
def getProcessedValues(varid,inputdir,inputfiles=None,dates=None,concataxis=0):

    """Get $varid values from processed data directory $inputdir."""

    if inputfiles is None:
        inputfiles = getInputfiles(varid,inputdir,inputfiles,dates)
    if len(inputfiles) == 0:
        print("Processed files: no acceptable input file for %s."%varid)
        return
    values_list = []
    if len(inputfiles) == 0:
        return np.array([])
    for file in inputfiles:
        fh = Dataset(file,'r')
        values_list.append(fh.variables[varid][:])
        fh.close()
    if len(values_list) == 0:
    	print("%s not found in processed files"%varid)
    	return None
    print("%s found in processed files"%varid)
    values_array = np.concatenate(values_list,axis=concataxis)
    return values_array

## Get values from original hourly output
def getSimulationValues(varid,inputdir,dates=None,subset='tropics',handle=None):

	"""Get $varid values from CAM/SPCAM history files in $inputdir.
	$dates is a tuple of time boundaries (YYYYMMDDhhmm <= inputfiles < YYYYMMDDhhmm)"""

	## Find valid inputfiles
	# Define search pattern based on file handle:
	pattern = "*.????-??-??-?????.nc"
	if handle is not None:
		pattern = "*.%s%s"%(handle,pattern[1:])
	# Search inputfiles
	inputfiles = []
	if dates is not None:
		dt_bnds = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8]),int(d[8:10]),int(d[10:12]))\
		for d in dates]
	for file in glob.glob(os.path.join(inputdir,pattern)):
		filename = os.path.basename(file)
		if dates is not None: 
			dt_info = filename.split('.')[-2].replace('-','')
			dt = datetime(int(dt_info[:4]),int(dt_info[4:6]),int(dt_info[6:8]))+\
				timedelta(seconds=int(dt_info[8:13]))
			if dt < dt_bnds[0] or dt >= dt_bnds[1]:
				continue
		inputfiles.append(file)
	inputfiles.sort()
	## Extract data

	# Abort if no inputfile match date range
	if len(inputfiles) == 0:
		print("    History files: no acceptable input file for %s.")%varid
		return

	values_list = []
	for file in inputfiles:
		fh = Dataset(file,'r')
		if varid in fh.variables.keys():
			values_list.append(fh.variables[varid][:])
	fh.close()

	# Abort if no inputfile contains varid
	if len(values_list) == 0:
		print("%s not found in history files"%varid)
		return 
	# If found
	print("%s found in history files"%varid)
	print("  Importing %s\n  from\n  %s\n  to\n  %s"%(varid,inputfiles[0],inputfiles[-1]))
	values = np.concatenate(values_list,axis=0)
	## Reduce domain
	if subset is not 'global':
		fh = Dataset(inputfiles[0],'r')
		dims = fh.variables[varid].dimensions
		ndims = len(dims)
		fh.close()
	if subset == 'tropics':
		indlat = np.argmax(np.array(dims) == 'lat')
		nlat = values.shape[indlat]
		lat_slice = slice(nlat//3,2*(nlat//3))
		lat_slice_multidim = [slice(None)]*(indlat)+[lat_slice]+[slice(None)]*(ndims-indlat-1)
		values = values[lat_slice_multidim]
	return values


def getValues(varid,dataroot,compset,subset,experiment,time_stride,time_type='A',dates=None):

	"""Wrapper around the other 'get*Values' functions that loads the data 
	for a given time resolution (stride), type and date range, for a compset
	and experiment.
	Primarily load data from CAM history files; if the variable does not exist,
	look at the preprocessed variables."""

	# Load history files settings and pick corresponding input directory/historyfile index
	settings = getCAMHistoryFilesSettings()
	handle = None
	for k in settings.keys():
		v = settings[k]
		if v[0] == time_stride and v[1] == time_type:
			handle = k
			pass
	if handle is None:
		print("No history file has the correct time_stride or averaging type.\n Rewrite function.")
	# Load inputdirectories
	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results,\
	 inputdir_fx = getInputDirectories(dataroot,compset,experiment)
	# Try history files in case this is an original varid
	values = getSimulationValues(varid,inputdir,dates,subset,handle)
	if values is not None:
		return values

	# Now look at preprocessed variables
	inputdir_processed = locals()["inputdir_processed_%s"%time_stride]
	values = getProcessedValues(varid,inputdir_processed,inputfiles=None,dates=dates)

	return values


## Reads in file lev_fx_* with the required information
## Returns a function which takes a surface pressure value
## and returns a vector of pressure values
def getPressureCoordinateFunction(input_lev_file):

	"""Returns a function that derives the pressure coordinate values
	from pressure values at the surface sigma level."""

	fh = Dataset(input_lev_file,'r')
	hyam = fh.variables['hyam'][:]
	hybm = fh.variables['hybm'][:]
	P0 = fh.variables['P0'][:]
	fh.close()
	def computeP(ps):
		if type(ps).__module__ == np.__name__:
			pres = np.add(np.multiply.outer(P0*np.ones(ps.shape),hyam),\
				np.multiply.outer(ps,hybm)) # In hPa,mbar
			return pres
		else:
			return P0*hyam+ps*hybm	# In hPa,mbar
	return computeP

