#!/bin/bash
#SBATCH --job-name="analyzePointWiseScalingOmega500TsPs"
#SBATCH --output="../logs/analyzePointWiseScalingOmega500TsPs.%j.%N.out"
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH --constraint=haswell

module load python/3.6-anaconda-4.4

# Convert notebook to python script
./ipynb2py36.sh analyzePointWiseScalingOmega500TsPs.ipynb analyzePointWiseScalingOmega500TsPs

# Edit options in script
sed -i='' 's/^dates =.*/dates = "185005010000","185006010000"/g' analyzePointWiseScalingOmega500TsPs.py

# Run script
python -u analyzePointWiseScalingOmega500TsPs.py &> ../logs/analyzePointWiseScalingOmega500TsPs_`date +"%Y%m%d-%H%M"`.log
