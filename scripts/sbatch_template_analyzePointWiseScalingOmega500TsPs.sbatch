#!/bin/bash
#SBATCH --mail-user=bfildier@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --job-name="analyzePointWiseScalingOmega500TsPs"
#SBATCH --output="../logs/analyzePointWiseScalingOmega500TsPs.%j.%N.out"
#SBATCH --account=m1517
##SBATCH --partition=regular
##SBATCH --partition=shared
#SBATCH --partition=debug
##SBATCH --qos=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
##SBATCH --cpus-per-task=30
#SBATCH --time=00:30:00
##SBATCH --time=9:00:00
#SBATCH --constraint=haswell
##SBATCH --constraint=knl,quad,cache

scriptname=
module load python/3.6-anaconda-4.4
python -u ${scriptname}.py > logs/${scriptname}.log

