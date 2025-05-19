#!/bin/bash
#SBATCH -N 5               # request one node
#SBATCH -n 100
#SBATCH -t 10:00:00	        # request 10 hours
#SBATCH -p workq          # in single partition (queue)
#SBATCH -A hpc_sel_001
#SBATCH -o 2023_damsel-%j.out-%N # optional, name of the stdout, using the job number (%j) and the hostname of the node (%N)
#SBATCH -e 2023_damsel-%j.err-%N # optional, name of the stderr, using job and hostname values
# below are job commands

date

module purge
srun -N5 -n100 singularity run -B /work,/project /project/containers/images/cms.08212024.sif cms 2023_damsel
# Mark the time it finishes.
date
exit 0
