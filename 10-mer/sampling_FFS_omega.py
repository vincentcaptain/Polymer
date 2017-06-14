from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from lammps import lammps
from mpi4py import MPI
from pylab import *
import random
import sys

"""
This script is to sample the initial condition for transition probability calculation. Input file has the
first argument for the starting point, and second argument is the size of sampling pool.
"""
pol = lammps(cmdargs = ["-sc", "none"])
pol.file("in.ljwall")
start = eval(sys.argv[1])
size = eval(sys.argv[2])
monomer = eval(sys.argv[3])
init_x = []
init_v = []
init_t = []
i = 0
pol.command("run 0 pre no post no")
com_current = pol.extract_compute("com", 0, 1)[0]
left = com_current < start
while i < size:
	com_current = pol.extract_compute("com", 0, 1)[0]
	if com_current > start and left:
		x_current = pol.gather_atoms("x", 1, 3)
		v_current = pol.gather_atoms("v", 1, 3)
		init_x.append(x_current[:3*monomer])
		init_v.append(v_current[:3*monomer])
		init_t.append(pol.extract_global("ntimestep", 0))
		i += 1
		left = not left
	if com_current < start and not left:
		left = not left
	pol.command("run 10")

savetxt("init_x.txt", init_x)
savetxt("init_v.txt", init_v)
savetxt("init_t.txt", init_t)
sys.exit()