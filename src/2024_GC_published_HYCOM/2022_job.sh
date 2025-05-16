#!/bin/bash
#SBATCH -N 4               # request one node
#SBATCH -n 80
#SBATCH -t 10:00:00	        # request 10 hours
#SBATCH -p workq          # in single partition (queue)
#SBATCH -A hpc_sel_001
#SBATCH -o 2022-%j.out-%N # optional, name of the stdout, using the job number (%j) and the hostname of the node (%N)
#SBATCH -e 2022-%j.err-%N # optional, name of the stderr, using job and hostname values
# below are job commands

date

module purge
srun -N4 -n80 singularity run -B /work,/project /project/containers/images/cms.08212024.sif cms 2022
# Mark the time it finishes.
date
exit 0
