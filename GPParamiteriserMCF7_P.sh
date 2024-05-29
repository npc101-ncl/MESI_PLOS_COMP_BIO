#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
#SBATCH --mem-per-cpu=5G
#

python GPParamiteriser.py slurm dataName:GPMCF7 data:zr75_data_MCF7only.xlsx twoStage:yes parallelAntStr:antString2MP_MCF7.txt methP:method:particle_swarm,swarm_size:100,iteration_limit:5000 stabAlow:200.0 
