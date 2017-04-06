from __future__ import division
from pylab import *
from lammps import lammps, PyLammps
from mpi4py import MPI
pol = lammps()
pol.file("in.ljwall")
polymer = PyLammps(ptr = pol)

# sampling range
sampling = [i / 10 + 27 for i in range(55)]

# Initialize
polymer.fix(3, "polymer spring tether 1 27 NULL NULL 0")
polymer.run(20000)
polymer.command("unfix 3")

def FFS(steps = 50):
	ptotal = 1
	stored_traj = []
	for x in sampling:
		Pass = 0
		Fail = 0
		while Pass < steps:
			
