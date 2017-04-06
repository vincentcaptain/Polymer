from __future__ import division
from lammps import PyLammps
from lammps import lammps
import sys
from mpi4py import MPI
comm = MPI.COMM_WORLD
print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())
p = lammps()
p.file("in.20mer")
i = 27
p.command("fix 3 polymer spring tether 1 33 NULL NULL 0")
n = p.extract_compute("com",0, 1)[0]
while n < 32:
	p.command("run 50000")
	n = p.extract_compute("com", 0, 1)[0]
	i += 50000
j = 0

p.command("unfix 3")
p.command("run 20000")
p.command("fix 3 polymer spring tether 1 27 NULL NULL 0")
while n > 28:
       p.command("run 50000")
       n = p.extract_compute("com", 0, 1)[0]
       j += 50000
print(i, j)