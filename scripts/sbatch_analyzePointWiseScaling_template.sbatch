#!/bin/bash
#SBATCH --mail-user=bfildier@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --job-name="analyzePointWiseScaling"
#SBATCH --output="../logs/analyzePointWiseScaling.%j.%N.out"
#SBATCH --account=m1517
#SBATCH --partition=debug
##SBATCH --qos=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --time=00:30:00
#SBATCH --constraint=haswell

scriptname=
module load python/3.6-anaconda-4.4
python -u ${scriptname}.py > logs/${scriptname}.log

