"""Module CAMsettings

Contains functions related to history files settings, time slices, etc.
"""

#---- modules ----#

import socket, os

from environmentAndDirectories import *

#---- Functions ----#

## Get time stride, type (A or I) and number or time slices for CAM history files
def getCAMHistoryFilesSettings():

	"""Returns a dict where keys are file handles h0, h1, ... and values are
	tuples (time stride, average type, number of time slices per file)."""

	historyFilesSettings = dict()
	historyFilesSettings['h0'] = ('1hr','A',1)
	historyFilesSettings['h1'] = ('day','A',1)

	return historyFilesSettings

## Returns a conversion factor to turn into hours
def inHours(resolution):
    
    """Returns the number of hours for a given time stride/resolution"""
    
    resolutions = ['day','1d','3hr','1hr']
    nhours = [24,24,3,1]
    conversion = dict(zip(resolutions,nhours))
    
    if resolution not in conversion.keys():
        print("Time resolution %s not known."%resolution)
        return
    
    return conversion[resolution]

## From a target and a source time resolutions, get the ratio
def timeResolutionRatio(target,source):
    
    """Ratio is target/source.
    If ratio < 1, send a warning and return None.
    If ratio > 1 but not integer, send a warning and return None."""
    
    n_t = inHours(target)
    n_s = inHours(source)
    
    ratio = n_t/n_s
    if ratio < 1:
        print("Required time resolution is too small.")
        return
    if int(ratio) != ratio:
        print("Required time resolution is not a multiple of the available resolution %s."%source)
        return
    
    return int(ratio)
    
## Get file handles from a list of CAM history file names
def handleCAMHistoryFiles(file):
    
    """Returns file handles h?."""
    
    getHandle = lambda f: f.split('/')[-1].split('.')[2]
    
    # If it's a single filename
    if file.__class__ is str:
        return getHandle(file)
    # If it's a list of filenames
    elif file.__class__ is list:
        handles = []
        for f in file:
            handles.append(f.split('/')[-1].split('.')[2])
        return handles
    
## Check whether the handle correspond to a file that can be
## used with the time resolution required
def isValidHandle(h,dt_target,settings):

    (dt_h,avg_h,n_t) = settings[h]
    dt_ratio = timeResolutionRatio(dt_target,dt_h)

    return dt_ratio is not None
