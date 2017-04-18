from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from lammps import lammps
from mpi4py import MPI
from pylab import *
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
	time = 0
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
		if com_current < target - 2:
			fail += 1
			r = random.randint(0, len(x) - 1)
			x_current[:3*monomers] = x[r]
			v_current[:3*monomers] = v[r]
			pol.scatter_atoms("x", 1, 3, x_current)
			pol.scatter_atoms("v", 1, 3, v_current)
		if Pass >= steps:
			break
		pol.command("run 1 pre no post no")
		time += 0.01
		print(time, Pass)
	ptotal = Pass/(Pass+fail)
	flux = Pass/time
	return init_x, init_v, ptotal, time


def FFS_cont(init_x, init_v, p_init, pos_init, target_series, monomers = 10, steps = 50):
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
				print(Pass)
			if com_current < pos_init:
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

def flux(pos, target):
	p = lammps(cmdargs = ["-sc", "none"])
	p.file("in.ljwall")
	p.command("fix 3 polymer spring tether 1 %d NULL NULL 0" % pos)
	p.command("run 5000")
	p.command("unfix 3")
	com = p.extract_compute("com", 0, 1)[0]
	if com > pos:
		cross = 1
		left = False
	else:
		cross = 0
		left = True
	t = 0
	while cross < target:
		p.command("run 50")
		com = p.extract_compute("com", 0, 1)[0]
		t += 0.5
		if (com > pos and left) or (com < pos and not left):
			cross += 1
			left = not left
		print(cross, com[1])
	return cross/time



# start = 27
# end = 33
# interval = 0.1
# monomer = 10
# steps = 50
# size = 900
# pol = lammps(cmdargs = ["-sc", "none"])
# pol.file("in.ljwall")
# # sampling range
# length = int((end - start)/interval)
# sampling = arange(start, end, 0.1)
# # Initialize
# pol.command("fix 3 polymer spring tether 1 26 NULL NULL 0")
# pol.command("run 5000 pre no post no")
# ptotal = 1
# x_init = []
# v_init = []
# com_y = []
# com_z = []
# while len(x_init) < size:
# 	pol.command("run 100 pre no post no")
# 	x_original = pol.gather_atoms("x", 1, 3)
# 	v_original = pol.gather_atoms("v", 1, 3)
# 	com_current = pol.extract_compute("com", 0, 1)
# 	if com_current[0] > 26:
# 		x_init.append(x_original[:3*monomer])
# 		v_init.append(v_original[:3*monomer])
# 		com_y.append(com_current[1])
# 		com_z.append(com_current[2])
# 	print(len(x_init))


# pol.command("unfix 3")
# pol.command("run 0")
# Q0 = FFS_init(x_init, v_init, start, monomer, size)
#Q1 = FFS_cont(Q0[0], Q0[1], Q0[2], start, sampling, monomer, size)
#np.savetxt("FFS_prob_and_com.txt", Q1)
#Q1 = FFS_cont(Q0[0], Q0[1], Q0[2], sampling, monomer, steps)
#np.savetxt("FFS_prob_and_com.txt", np.r_[sampling, Q1])
#sys.exit()


