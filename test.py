from __future__ import division
from lammps import PyLammps
from lammps import lammps
p = lammps()
p.file("in.ljwall")
i = 0
while i < 16:
	j = i / 2
	p.command("fix 3 polymer spring tether 5 %d NULL NULL 0" % j)
	p.command("run 500")
	p.command("unfix 3")
	i += 1
