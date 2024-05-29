#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
#SBATCH --mem-per-cpu=5G
#

python -u GPParamiteriser.py slurm scaleIn:scales.csv configFile:scalemodelkeys.json dataName:GPMCF7 copys:20 paramTrans:parameters.csv paramOut:parameters_ZR75.csv scaleOut:scales_ZR75.csv cmdRec:comandRecordZR75.p TCOutputs:PETimeCourses_ZR75.csv DBOutputs:DBTimeCourses_ZR75.csv elemOut:elements_ZR75.p data:zr75_data_ZR75only.xlsx twoStage:yes parallelAntStr:antString2MP_ZR75.txt methP:method:particle_swarm,swarm_size:100,iteration_limit:5000 stabAlow:500.0 parameterUB:100000.0 parameterLB:0.00001 
