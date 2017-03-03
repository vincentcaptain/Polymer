import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpi4py import MPI
me = MPI.COMM_WORLD.Get_size()
from lammps import lammps, PyLammps
lmp1 = lammps()
lmp2 = lammps()
Lattice = PyLammps(ptr = lmp1)
polymer = PyLammps(ptr=lmp2)

# Eliminate repeated points and convert them to normal coordinate
def position_data_modify(x, y, z):
	i = 0
	j = 1
	n = 0
	while i < len(x):
	   while j < len(x):
	      if x[i] == x[j]:
	         if y[i] == y[j] and z[i] == z[j]:
	            x = x[:j]+x[j+1:]
	            y = y[:j]+y[j+1:]
	            z = z[:j]+z[j+1:]
	            j = j - 1
	      j += 1
	   i += 1
	   j = i+1

def CheckRepeats(x, y, z):
	i = 0
	j = 1
	n = 0
	while i < len(x):
	   while j < len(x):
	      if x[i] == x[j]:
	         if y[i] == y[j] and z[i] == z[j]:
	            return False
	      j += 1
	   i += 1
	   j = i+1
	return True