#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
#SBATCH --mem-per-cpu=5G
#

python GPParamiteriser.py slurm dataName:redPaper id:zr75 twoStage:yes antStr:antStringRedM.txt preAntStr:antStringRedF.txt parallelAntStr:antStringRedMP.txt stabAlow:200.0 methP:method:particle_swarm,swarm_size:100,iteration_limit:5000 scaleRefs:PEScaleRefs_zr75.csv paramOut:parameters_zr75.csv scaleOut:scales_zr75.csv TCOutputs:PETimeCourses_zr75.csv DBOutputs:DBTimeCourses_zr75.csv
