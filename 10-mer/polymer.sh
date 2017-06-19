#!/bin/bash
# Job name:
#SBATCH --job-name=sampling
#
# Account:
#SBATCH --account=co_noneq
#
# Partition:
#SBATCH --partition=savio2
#
# Request one node:
#SBATCH --nodes=1
#
# Processors per task:
#SBATCH --cpus-per-task=4
#
# Wall clock limit:
#SBATCH --time=80:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=vincentcaptain@berkeley.edu
## Command(s) to run:
mpirun -np 4 python2 sampling_FFS_omega.py 25.4 900 10
