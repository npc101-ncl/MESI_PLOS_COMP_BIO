#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
#SBATCH --mem-per-cpu=5G
#

python GPParamiteriser.py slurm dataName:paper twoStage:yes parallelAntStr:antString2MP.txt methP:method:particle_swarm,swarm_size:100,iteration_limit:5000 stabAlow:200.0 ratioWeight:2
