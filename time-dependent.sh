#!/bin/bash
# Job name:
#SBATCH --job-name=ljwall_24_34.5
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
#SBATCH --cpus-per-task=20
#
# Wall clock limit:
#SBATCH --time=40:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=vincentcaptain@berkeley.edu
## Command(s) to run:
mpirun -np 20 python2 test.py
