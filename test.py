from __future__ import division
from lammps import PyLammps
from lammps import lammps
import sys
from mpi4py import MPI
comm = MPI.COMM_WORLD
print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())
p = lammps()
p.file("in.ljwall")
i = 0
while i < 2:
	p.command("fix 3 polymer spring tether 5 %d NULL NULL 0" % i)
	p.command("run 1000000")
	p.command("unfix 3")
	i += 2
sys.exit()
