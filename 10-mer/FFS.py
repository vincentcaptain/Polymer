from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from lammps import lammps
from mpi4py import MPI
import random
import sys

# def main(argv):
# 	if len(argv) != 0 and len(argv) != 5:
# 		print "input should only include starting point, end point, interval, monomer number and sampling steps"
# 		sys.exit()
# 	if len(argv) == 5:
# 		start = eval(sys.argv[0])
# 		end = eval(sys.argv[1])
# 		interval = eval(sys.argv[2])
# 		monomer = eval(sys.argv[3])
# 		steps = eval(sys.argv[4])
# 	else:
def FFS_init(x, v, target, monomers = 10, steps = 50):
	init_v = []
	init_x = []
	Pass = 0
	fail = 0
	while Pass < steps:
		# Gather all atom position
		x_current = pol.gather_atoms("x", 1, 3)
		v_current = pol.gather_atoms("v", 1, 3)
		com_current = pol.extract_compute("com", 0, 1)[0]
		if com_current > target:
			Pass += 1
			init_v.append(v_current[:3*monomers])
			r = random.randint(0, len(x) - 1)
			init_x.append(x_current[:3*monomers])
			x_current[:3*monomers] = x[r]
			v_current[:3*monomers] = v[r]
			pol.scatter_atoms("x", 1, 3, x_current)
			pol.scatter_atoms("v", 1, 3, v_current)
			print(Pass)
		if com_current < target - 1:
			fail += 1
			r = random.randint(0, len(x) - 1)
			x_current[:3*monomers] = x[r]
			v_current[:3*monomers] = v[r]
			pol.scatter_atoms("x", 1, 3, x_current)
			pol.scatter_atoms("v", 1, 3, v_current)
		if Pass >= steps:
			break
		pol.command("run 1 pre no post no")
	ptotal = Pass/(Pass+fail)
	return init_x, init_v, ptotal


def FFS_cont(init_x, init_v, p_init, target_series, monomers = 10, steps = 50):
	p_collection = []
	for item in target_series:
		Pass = 0
		fail = 0
		cont_x = []
		cont_v = []
		while Pass < steps:
			x_current = pol.gather_atoms("x", 1, 3)
			v_current = pol.gather_atoms("v", 1, 3)
			com_current = pol.extract_compute("com", 0, 1)[0]
			if com_current > item:
				Pass += 1
				cont_v.append(v_current[:3*monomers])
				cont_x.append(x_current[:3*monomers])
				r = random.randint(0, len(init_x) - 1)
				x_original[:3*monomers] = init_x[r]
				v_original[:3*monomers] = init_v[r]
				pol.scatter_atoms("x", 1, 3, x_original)
				pol.scatter_atoms("v", 1, 3, v_original)
			if com_current < item - 1:
				fail += 1
				r = random.randint(0, len(init_x) - 1)
				x_original[:3*monomers] = init_x[r]
				v_original[:3*monomers] = init_v[r]
				pol.scatter_atoms("x", 1, 3, x_original)
				pol.scatter_atoms("v", 1, 3, v_original)
			if Pass >= steps:
				break
			pol.command("run 1 pre no post no")
		p_init *= Pass/(Pass+fail)
		p_collection.append(Pass/(Pass+fail))
		init_x = cont_x
		init_v = cont_v
	return p_init, p_collection


start = 27
end = 32
interval = 0.1
monomer = 10
steps = 10
pol = lammps()
pol.file("in.ljwall")
# sampling range
length = int((end - start)/interval)
sampling = [i*interval + start for i in range(end)]
# Initialize
pol.command("fix 3 polymer spring tether 1 %s NULL NULL 0" % start)
pol.command("run 5000 pre no post no")
ptotal = 1
x_init = []
v_init = []
while len(x_init) < 50:
	pol.command("run 20 pre no post no")
	x_original = pol.gather_atoms("x", 1, 3)
	v_original = pol.gather_atoms("v", 1, 3)
	com_current = pol.extract_compute("com", 0, 1)[0]
	if com_current > 25.5:
		x_init.append(x_original[:3*monomer])
		v_init.append(v_original[:3*monomer])

pol.command("unfix 3")
pol.command("run 0")
Q0 = FFS_init(x_init, v_init, start, monomer, steps)
#Q1 = FFS_cont(Q0[0], Q0[1], Q0[2], sampling, monomer, steps)
#np.savetxt("FFS_prob_and_com.txt", np.r_[sampling, Q1])



