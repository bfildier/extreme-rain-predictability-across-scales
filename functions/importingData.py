"""Module importingData

Contains
- functions to import variables between specific dates at various time
resolutions, from model output files or processed files.
If the variable exists in CAM, directly import it from CAM history files. If
not, import it from the processed variables.
- functions to import fixed quantities (pressure coordinate function, areacella)
"""

#---- Modules ----#

import glob
import sys,os
import string
import operator
import numpy as np
import dask.array as da
from datetime import date, datetime, timedelta
from netCDF4 import Dataset

## Own functions
# currentpath = os.path.dirname(os.path.realpath(__file__))

from environmentAndDirectories import *
from CAMsettings import *
from slicingAndSubsetting import *

#---- Functions ----#

## inputfiles is a list of inputfiles to open simultaneously.
## dates is a pair of date strings in the form 'YYYYMMDD'
## - if both arguments 'inputfiles' and 'dates' are NA, 
##   use all files $inputdir/$varid_*
## - if inputfiles is given, use this list in priority
## - if inputfiles is not given but dates is given, use corresponding files
## between those dates
## Get processed file names in correct format 
def getProcessedFiles(varid,inputdir,inputfiles=None,dates=None):

	"""Based on a variable ID and input directory, output a list of corresponding 
	filenames for preprocessed files.
	 - $inputfiles are in format $varid_*_date1-date2.nc
	 (add directory to filenames if not present)
	 - dates is a pair of dates in format YYYYmmddHHMM
	 (YYYYmmddHHMM <= accepted files < YYYYmmddHHMM)"""

	# print("enter getProcessedFiles")

	## Get the list of input files
	if inputfiles is None:
		if dates is None:    # Get all matching files
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
		if dates is None:
			dates=[]
		print("No processed file found for varid %s and dates %s in %s."%\
			(varid,'-'.join(dates),inputdir))
	return inputfiles

## Get values from processed data files ($dataroot/preprocessed/$case/$freq/*)
def getProcessedValues(varid,inputdir,inputfiles=None,dates=None,concataxis=0,daskarray=True):

    """Get $varid values from processed data directory $inputdir."""

    cn = da if daskarray else np
    
    if inputfiles is None:
        inputfiles = getProcessedFiles(varid,inputdir,inputfiles,dates)
    if len(inputfiles) == 0:
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
    # print("%s found in processed files"%varid)
    print("Importing %s from %d processed files between %s and %s"%
		(varid,len(inputfiles),inputfiles[0].split('.')[-2],
			inputfiles[-1].split('.')[-2]))
    values_array = cn.concatenate(values_list,axis=concataxis)
    return values_array

def getSimulationFiles(varid,inputdir,inputfiles=None,dates=None,handle=None):

	"""Based on a variable ID and input directory, return a list of corresponding
	filenames for simulation files (CAM history files).
	- $inputfiles are in format *.h?.YYYY-mm-dd-SSSSS.nc
	 (add directory to filenames if not present)
	- dates is a pair of dates in format YYYYmmddHHMM
	(YYYYmmddHHMM <= accepted files < YYYYmmddHHMM)"""

	# print("enter getSimulationFiles")

	# Define search pattern based on file handle:
	pattern = "*.????-??-??-?????.nc"
	if handle is not None:
		pattern = "*.%s%s"%(handle,pattern[1:])
	# Search inputfiles
	if inputfiles is None:
		inputfiles = []
		if dates is not None:
			dt_bnds = [datetime(int(d[:4]),int(d[4:6]),int(d[6:8]),int(d[8:10]),int(d[10:12]))\
			for d in dates]
		for file in glob.glob(os.path.join(inputdir,pattern)):
			if dates is not None: 
				dt_info = file.split('.')[-2].replace('-','')
				dt = datetime(int(dt_info[:4]),int(dt_info[4:6]),int(dt_info[6:8]))+\
					timedelta(seconds=int(dt_info[8:13]))
				if dt < dt_bnds[0] or dt >= dt_bnds[1]:
					continue
			inputfiles.append(file)
	else:
		## Append dirname to all files if necessary
		inputfiles = [os.path.join(inputdir,f) for f in inputfiles]
	
	inputfiles.sort()
	if len(inputfiles) == 0:
		if dates is None:
			dates=['undefined']
		print("No simulation history file found for dates %s in %s."%\
			('-'.join(dates),inputdir))
	return inputfiles

## Get values from original hourly output
def getSimulationValues(varid,inputdir,dt='day',inputfiles=None,dates=None,subsetname='tropics',handle=None,daskarray=True):

	"""Get $varid values from CAM/SPCAM history files in $inputdir.
	dates is a pair of dates in format YYYYmmddHHMM
	(YYYYmmddHHMM <= accepted files < YYYYmmddHHMM).
	Assumes the first dimension is time."""
	
	cn = da if daskarray else np

	# Find valid inputfiles
	if inputfiles is None:
		inputfiles = getSimulationFiles(varid,inputdir,inputfiles,dates,handle)
	# Abort if no inputfile matches date range
	if len(inputfiles) == 0:
		return
    
	# Check time resolution compatibility
	if handle is None:
		handle = handleCAMHistoryFiles(inputfiles[0])
	settings = getCAMHistoryFilesSettings()
	if not isValidHandle(handle,dt,settings):
		return 
    
	# Import data
	dt_ratio = timeResolutionRatio(dt,settings[handle][0])
    
	vals_within_dt = []
	values_list = []
	for file in inputfiles:
		fh = Dataset(file,'r')
		if varid in fh.variables.keys():
			vals_within_dt.append(fh.variables[varid][:])
		if len(vals_within_dt) == dt_ratio:
			values_list.append(cn.mean(cn.concatenate(vals_within_dt,axis=0),
									   axis=0,
									   keepdims=True))
			vals_within_dt = []
		fh.close()

	# Abort if no inputfile contains varid
	if len(values_list) == 0:
		print("%s not found in history files"%varid)
		return 
    
	# If found
	# print("%s found in history files"%varid)
	# print("Importing %s from %d files\n  from\n  %s\n  to\n  %s"%
	# 	(varid,len(inputfiles),inputfiles[0],inputfiles[-1]))
	print("Importing %s from %d history files between %s and %s"%
		(varid,len(inputfiles),inputfiles[0].split('.')[-2],
			inputfiles[-1].split('.')[-2]))
	values = cn.concatenate(values_list,axis=0)
	values = reduceDomain(values,subsetname,inputfiles[0],varid)

	return values

# Main function to extract data values from processed files or history files
def getValues(varid,compset,subsetname,experiment,time_stride,time_type='A',dates=None,handle=None,daskarray=True):

	"""Wrapper around the other 'get*Values' functions that loads the data 
	for a given time resolution (stride), type and date range, for a compset
	and experiment.
	Primarily load data from CAM history files; if the variable does not exist,
	look at the preprocessed variables."""

	# print("enter getValues")

	# Load inputdirectories
	inputdir, inputdir_processed_day, inputdir_processed_1hr, inputdir_results,\
	 inputdir_fx = getInputDirectories(compset,experiment)
	
	# Try history files in case this is an original varid
	values = getSimulationValues(varid,inputdir,time_stride,None,dates,subsetname,handle,daskarray=daskarray)
	if values is not None:
		return values

	# Now look at preprocessed variables
	inputdir_processed = locals()["inputdir_processed_%s"%time_stride]
	values = getProcessedValues(varid,inputdir_processed,inputfiles=None,dates=dates,daskarray=daskarray)

	return values

## Obtains the 'computeP' function for a given simulation
def getPressureCoordinateFunction(input_lev_file,levdim=1):

	"""Reads in file lev_fx_* with the required information (hyam, hybm, P0).
	Returns a function that derives the pressure coordinate values on sigma levels
	from arrays of surface pressure values."""

	fh = Dataset(input_lev_file,'r')
	hyam = fh.variables['hyam'][:]
	hybm = fh.variables['hybm'][:]
	P0 = fh.variables['P0'][:]
	fh.close()
	def computeP(ps):
		if ps.__class__ == np.ndarray:
			pres = np.add(np.multiply.outer(P0*np.ones(ps.shape),hyam),\
				np.multiply.outer(ps,hybm)) # In hPa,mbar
			return np.moveaxis(pres,-1,levdim)
		elif ps.__class__ == da.core.Array:
			chunks = ps.chunks
			newchunks = chunks[:levdim]+((hyam.size,),)+chunks[levdim:]
			# pres = da.map_blocks(computeP,ps,dtype=float).compute() ## crashes
			pres = computeP(ps.compute())
			return(da.from_array(pres,chunks=newchunks))
		else:
			return P0*hyam+ps*hybm	# In hPa,mbar
	return computeP

