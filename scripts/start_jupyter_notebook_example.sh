#!/bin/bash
#SBATCH --job-name="jupyter-notebook"
#SBATCH --output="jupyter-notebook.%j.%N.out"
#SBATCH --partition=debug
#SBATCH --nodes=32
##SBATCH --nodes=32
##SBATCH --ntasks-per-node=64
#SBATCH --export=ALL
#SBATCH --time=00:20:00
#SBATCH --constraint=haswell


module load python/3.6-anaconda-4.4

jupyter notebook --no-browser --ip=* --port=8889
