#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
#SBATCH --mem-per-cpu=5G
#

python GPParamiteriser.py slurm configFile:toyScriptSub.json