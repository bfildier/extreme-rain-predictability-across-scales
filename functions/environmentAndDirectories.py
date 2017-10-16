"""Module environmentAndDirectories

Contains functions to define directory structure, paths and other settings.
"""

#---- modules ----#

import socket, os

#---- Own functions ----#
# currentpath = os.path.dirname(os.path.realpath(__file__))
## Host name
hostname = socket.gethostname()
if hostname == "jollyjumper":
	## Parent directory to load data
	dataroot = "/Users/bfildier/Data"
elif "edison" in hostname or "cori" in hostname:
	dataroot = "/global/cscratch1/sd/bfildier/lawrencium_runs"
# sys.path.insert(0,os.path.join(os.path.dirname(currentpath),'functions'))

#---- Functions ----#

## Define input directories for different platforms
def getInputDirectories(compset,experiment):

	"""Arguments:
		- compset is FSPCAMm_AMIP or FAMIC5
		- experiment is piControl of abrupt4xCO2
	Returns inputdir, inputdir_processed_day, inputdir_processed_1hr,
		inputdir_results, inputdir_fx."""

	case = "bf_%s_%s"%(compset,experiment)
	if hostname == "jollyjumper":
		inputdir = os.path.join(dataroot,"simulations",case)
		inputdir_processed_day = os.path.join(dataroot,'preprocessed',case,'day')
		inputdir_processed_1hr = os.path.join(dataroot,'preprocessed',case,'1hr')
	elif "edison" in hostname or "cori" in hostname:
		inputdir = os.path.join(dataroot,'archive',case,"atm/hist")
		inputdir_processed_day = os.path.join(dataroot,'preprocessed',
			case,'day')
		inputdir_processed_1hr = os.path.join(dataroot,'preprocessed',
			case,'1hr')
	inputdir_results = os.path.join(dataroot,'results')
	inputdir_fx = os.path.join(dataroot,'preprocessed/allExperiments/fx')

	return inputdir, inputdir_processed_day, inputdir_processed_1hr,\
		inputdir_results, inputdir_fx

	