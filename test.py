from __future__ import division
from lammps import PyLammps
from lammps import lammps
import sys
p = lammps()
p.file("in.ljwall")
i = 46
while i < 60:
	p.command("fix 3 polymer spring tether 5 %d NULL NULL 0" % i)
	p.command("run 500")
	p.command("unfix 3")
	i += 2
sys.exit()
