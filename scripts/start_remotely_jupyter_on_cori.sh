#!/bin/bash

remote_port=$1
local_port=$2

#-----------------#
#	Directories	  #
#-----------------#

LOCAL_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# REMOTE_WORKDIR="/global/cscratch1/sd/bfildier/dataAnalysis/"\
# "extreme-rain-predictability-across-scales" # Doesn't work
REMOTE_WORKDIR="/global/homes/b/bfildier/code/dataAnalysis/"\
"extreme-rain-predictability-across-scales/"
# REMOTE_WORKDIR="/global/homes/b/bfildier/code/jupyter_dask_examples/"
REMOTE_SCRIPT_DIR="/global/homes/b/bfildier/scripts/jupyter/"
SCRIPT_NAME="start_jupyter_notebook.sh"

#------------------------#
#	Run specifications   #
#------------------------#

machine="cori"
constraint="haswell"
partition="regular"
nodes=2
time=04:00:00

#---------------------#
#	Launch notebook   #
#---------------------#

ssh -Y ${machine}.nersc.gov <<EOF
# Edit script
echo "cd ${REMOTE_SCRIPT_DIR}"
cd ${REMOTE_SCRIPT_DIR}
sed -i "s/#SBATCH --partition=.*/#SBATCH --partition=${partition}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --nodes=.*/#SBATCH --nodes=${nodes}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --time=.*/#SBATCH --time=${time}/" ${SCRIPT_NAME}
sed -i '' "s/port=.*/port=${remote_port}/" ${SCRIPT_NAME}
# Start notebook
echo "cd ${REMOTE_WORKDIR}"
cd ${REMOTE_WORKDIR}
echo "Execute jupyter launch script"
touch jupyter.log
# ${REMOTE_SCRIPT_DIR}${SCRIPT_NAME} 2>&1 | tee jupyter.log && exit 0
${REMOTE_SCRIPT_DIR}${SCRIPT_NAME} &> >(tee jupyter.log)
EOF

scp ${machine}.nersc.gov:${REMOTE_WORKDIR}jupyter.log ${LOCAL_DIR}/
ssh ${machine}.nersc.gov "rm ${REMOTE_WORKDIR}jupyter.log"

remote_port=`cat ${LOCAL_DIR}/jupyter.log | grep "The Jupyter Notebook is running at" | cut -d: -f6 | cut -d\/ -f1` && echo "Press ^C"

ssh -Y -N -f -L localhost:${local_port}:localhost:${remote_port} ${machine}.nersc.gov

echo "Connecting local port ${local_port} to remote port ${remote_port}"

open http://localhost:${local_port}

exit 0



