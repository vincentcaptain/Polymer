#!/bin/bash
# Job name:
#SBATCH --job-name=ljwall
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
#SBATCH --time=20:00:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=vincentcaptain@berkeley.edu
## Command(s) to run:
python2 test.py
