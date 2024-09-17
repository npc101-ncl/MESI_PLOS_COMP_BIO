#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peter.clark@ncl.ac.uk
#

python scaleFinderXParam.py tor5 meth:method:particle_swarm,swarm_size:550,iteration_limit:10000 slurm isFullBGR mTorConst ICS globalChase
