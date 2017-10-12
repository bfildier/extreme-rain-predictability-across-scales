#!/bin/bash

#-----------------#
#	Directories	  #
#-----------------#

REMOTE_JUPYTER_DIR="/global/cscratch1/sd/bfildier/jupyter_test/"
REMOTE_SCRIPT_DIR="/global/homes/b/bfildier/scripts/jupyter/"
SCRIPT_NAME="start_jupyter_notebook.sh"

#------------------------#
#	Run specifications   #
#------------------------#

machine="cori"
constraint="haswell"
partition="debug"
nodes=32
time=00:20:00

#---------------------#
#	Launch notebook   #
#---------------------#

ssh -Y ${machine}.nersc.gov << EOF
# Edit script
cd ${REMOTE_SCRIPT_DIR}
sed -i "s/#SBATCH --partition=.*/#SBATCH --partition=${partition}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --nodes=.*/#SBATCH --nodes=${nodes}/" ${SCRIPT_NAME}
sed -i "s/#SBATCH --time=.*/#SBATCH --time=${time}/" ${SCRIPT_NAME}
# Start notebook
cd ${REMOTE_JUPYTER_DIR}
${REMOTE_SCRIPT_DIR}${SCRIPT_NAME}
EOF

ssh -Y -N -f -L localhost:8889:localhost:8889 cori.nersc.gov

open http://localhost:8889

exit 0



