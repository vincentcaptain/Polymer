from __future__ import division
from lammps import PyLammps
from lammps import lammps
from joblib import Parallel, delayed
import multiprocessing
i = 0
def run(i):
	p = lammps()
	p.file("in.ljwall")
	p.command("fix 3 polymer spring tether 5 %d NULL NULL 0" % i)
	p.command("run 50")
	#p.command("unfix 3")
	return p.extract_global("natoms", 1)


num_cores = multiprocessing.cpu_count()

results = Parallel(n_jobs = num_cores)(delayed(run)(i) for i in range(1, 10))