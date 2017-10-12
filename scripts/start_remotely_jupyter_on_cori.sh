#!/bin/bash

#-----------------#
#	Directories	  #
#-----------------#

# REMOTE_WORKDIR="/global/cscratch1/sd/bfildier/dataAnalysis/"\
# "extreme-rain-predictability-across-scales"
REMOTE_WORKDIR="/global/homes/b/bfildier/code/dataAnalysis/"\
"extreme-rain-predictability-across-scales/"
REMOTE_SCRIPT_DIR="/global/homes/b/bfildier/scripts/jupyter/"
SCRIPT_NAME="start_jupyter_notebook.sh"

#------------------------#
#	Run specifications   #
#------------------------#

machine="cori"
constraint="haswell"
partition="debug"
nodes=32
time=00:30:00

#---------------------#
#	Launch notebook   #
#---------------------#

ssh -Y ${machine}.nersc.gov << EOF
# Edit script
echo "cd ${REMOTE_SCRIPT_DIR}"
cd ${REMOTE_SCRIPT_DIR}
sed -i "s/#SBATCH --partition=.*/#SBATCH --partition=${partition}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --nodes=.*/#SBATCH --nodes=${nodes}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --time=.*/#SBATCH --time=${time}/" ${SCRIPT_NAME}
# Start notebook
echo "cd ${REMOTE_WORKDIR}"
cd ${REMOTE_WORKDIR}
echo "Execute jupyter launch script"
${REMOTE_SCRIPT_DIR}${SCRIPT_NAME}
EOF

ssh -Y -N -f -L localhost:8889:localhost:8889 cori.nersc.gov

open http://localhost:8889

exit 0



