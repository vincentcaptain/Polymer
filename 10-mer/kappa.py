from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from lammps import lammps
from mpi4py import MPI
from pylab import *
import random
import sys


start = 29.99
end = 30.01
interval = 0.1
monomer = 10
steps = 50
size = 64
pol = lammps(cmdargs = ["-sc", "none"])
pol.file("in.ljwall")
# sampling range

# Initialize
pol.command("fix 3 polymer spring tether 1 36 NULL NULL 0")
pol.command("run 1400")
ptotal = 1
x_init = []
v_init = []
com_y = []
com_z = []
pol.command("unfix 3")
pol.command("fix 3 polymer spring tether 5 30 NULL NULL 0")
while len(x_init) < size:
	pol.command("run 100")
	x_original = pol.gather_atoms("x", 1, 3)
	v_original = pol.gather_atoms("v", 1, 3)
	com_current = pol.extract_compute("com", 0, 1)
	if start <= com_current[0] <= end:
		x_init.append(x_original[:3*monomer])
		v_init.append(v_original[:3*monomer])
		com_y.append(com_current[1])
		com_z.append(com_current[2])
	print(len(x_init))


pol.command("unfix 3")


def kappa(left, right, left_end, right_end, x, v, target, monomers):
	left_total, right_total, left_recross, right_recross = 0, 0, 0, 0
	com_current = 0
	k = []
	while left_recross < target or right_recross < target:
		direction = push(left, right)
		if direction == "L":
			pol.command("run 3000")
			if check_left(left_end, right_end):
				left_total += 1
			else:
				right_total += 1
				right_recross += 1
		else:
			pol.command("run 3000")
			if check_left(left_end, right_end):
				left_total += 1
				left_recross += 1
			else:
				right_total += 1
		k.append(1-(left_recross+right_recross)/(left_total+right_total))
		r = random.randint(0, len(x) - 1)
		x_current = pol.gather_atoms("x", 1, 3)
		v_current = pol.gather_atoms("v", 1, 3)
		x_current[:3*monomers] = x[r]
		v_current[:3*monomers] = v[r]
		pol.scatter_atoms("x", 1, 3, x_current)
		pol.scatter_atoms("v", 1, 3, v_current)
		print(left_recross, right_recross)
	return k

def check_left(left, right):
	"""
	check if the com of polymer falls into left bin (T) or the right bin (F)
	"""
	com = pol.extract_compute("com", 0, 1)[0]
	while left <= com <= right:
		pol.command("run 300")
		com = pol.extract_compute("com", 0, 1)[0]
	if com < left:
		return True
	else:
		return False

	

def push(left, right):
	"""
	push the polymer to either side of the well
	"""
	com = pol.extract_compute("com", 0, 1)[0]
	while left <= com <= right:
		pol.command("run 1")
		com = pol.extract_compute("com", 0, 1)[0]
	if com < left:
		return "L"
	else:
		return "R"

result = kappa(start, end, 27.5, 32.5, x_init, v_init, size, monomer)
np.savetxt("kappa.txt", result)
