"""@package docstring
Module environmentAndDirectories

Contains functions to define directory structure and related history files 
settings.
"""

import socket, os

## Own functions
currentpath = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

###--- Functions ---###

# Define input directories for different platforms
def getInputDirectories(dataroot,compset,experiment):

	"""Define input directories for different platforms"""

	hostname = socket.gethostname()
	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,"simulations",case)
		inputdir_processed_day = os.path.join(dataroot,'preprocessed',case,'day')
		inputdir_processed_1hr = os.path.join(dataroot,'preprocessed',case,'1hr')
	elif "edison" in hostname or "cori" in hostname:
		inputdir = os.path.join(dataroot,'archive',case,"atm/hist")
		inputdir_processed_day = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'day')
		inputdir_processed_1hr = os.path.join(os.path.dirname(currentpath),'preprocessed',
			case,'1hr')
	inputdir_results = os.path.join(os.path.dirname(currentpath),'results')
	inputdir_fx = os.path.join(dataroot,'preprocessed/allExperiments/fx')

	return inputdir, inputdir_processed_day, inputdir_processed_1hr,\
		inputdir_results, inputdir_fx

def getCAMHistoryFilesSettings():

	"""Get time stride, type (A or I) and number or time slices per file."""

	historyFilesSettings = dict()
	historyFilesSettings['h0'] = ('1hr','A',1)
	historyFilesSettings['h1'] = ('day','A',1)

	return historyFilesSettings