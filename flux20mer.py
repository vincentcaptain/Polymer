from __future__ import division
from lammps import lammps
import sys
from mpi4py import MPI
p = lammps()
p.file("in.20mer")
i = 0
p.command("fix 3 polymer spring tether 1 30 NULL NULL 0")
n = p.extract_compute("com", 0, 1)[0]
j = 0
p.command("run 15000000")
sys.exit()
