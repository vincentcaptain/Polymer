from __future__ import division
from lammps import PyLammps
from lammps import lammps
import sys
from mpi4py import MPI
comm = MPI.COMM_WORLD
print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())
p = lammps()
p.file("in.ljwall")
i = 37
p.command("fix 3 polymer smoothforce 27 33 0.15 3")
while i > 22:
        p.command("fix 4 polymer spring tether 2 %d NULL NULL 0" % i)
        p.command("run 15000000")
        p.command("unfix 4")
	i -= 1
#while i < 60:
#	p.command("fix 3 polymer spring tether 3 %d NULL NULL 0" % i)
#        p.command("run 500000")
#        p.command("unfix 3")
#        i += 1
sys.exit()
