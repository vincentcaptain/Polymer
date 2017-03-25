from __future__ import division
from lammps import PyLammps
from lammps import lammps
import sys
from mpi4py import MPI
comm = MPI.COMM_WORLD
print "Proc %d out of %d procs" % (comm.Get_rank(),comm.Get_size())
p = lammps()
p.file("in.ljwall")
#i = 0
#while i < 24:
#	p.command("fix 3 polymer spring tether 3 %d NULL NULL 0" % i)
#	p.command("run 500000")
#	p.command("unfix 3")
#	i += 1
#i = 35
#while i < 60:
#	p.command("fix 3 polymer spring tether 3 %d NULL NULL 0" % i)
#        p.command("run 500000")
#        p.command("unfix 3")
#        i += 1
i = 24
while i < 35:
        p.command("fix 3 polymer spring tether 1 %d NULL NULL 0" % i)
        p.command("run 5000000")
        p.command("unfix 3")
        i += 0.5
sys.exit()
