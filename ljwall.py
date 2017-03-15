from __future__ import division
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random
# from mpi4py import MPI
# me = MPI.COMM_WORLD.Get_size()
from lammps import lammps, PyLammps
from joblib import Parallel, delayed
import multiprocessing
# l = []
# i = 90
# while i > 0:
# 	polymer.fix(3, "polymer spring tether", 10, 30, 15, 15, i/3)
# 	polymer.run(50000)
# 	l += [polymer.runs[90 - i][0][0][10000:]]
# 	polymer.unfix(3)
# 	i -= 1
# u = []
# j = 0
# while j < len(l):
# 	u += [[np.mean(l[j]), (90 - j)/3]]
# 	j += 1
# np.savetxt("all distribution", l)
# np.savetxt("mean", u)

# Get the mean force on it
def data(i):
	lmp1 = lammps()
	polymer = PyLammps(ptr=lmp1)
	x = 60
	y = 30
	z = 30
	t = 1
	polymer.units("lj")
	polymer.dimension(3)
	polymer.atom_style("bond")
	polymer.bond_style("harmonic")
	polymer.pair_style("lj/cut", 3)
	polymer.read_data("data.polymer")
	polymer.region("void cylinder x", 15, 15, 2, 29, 31)
	polymer.pair_coeff(1, 2, 2.5, 3)
	polymer.pair_coeff(1, 3, 2.5, 1.12)
	polymer.pair_coeff(2, 3, 2.5, 1.12)
	polymer.velocity("all create", t, 97287)
	polymer.group("polymer type", 1, 2)
	polymer.group("first type", 1)
	polymer.region("box block", 0, x, 0, y, 0, z)
	polymer.region("spherein sphere", 29, 15, 15, 2)
	polymer.region("boxin block", 27, 29, 13, 17, 13, 17)
	polymer.region("hemiin intersect", 2, "boxin spherein")
	x0 = polymer.atoms[0].position[0]
	y0 = polymer.atoms[0].position[1]
	z0 = polymer.atoms[0].position[2]
	r = lambda x0, y0, z0: np.sqrt((x0 - 29)**2 + (y0 - 15)**2 + (z0 - 15)**2)
	fx = lambda x0, y0, z0: 5*(x0-29)/r(x0, y0, z0)
	fy = lambda x0, y0, z0: 5*(y0-15)/r(x0, y0, z0)
	fz = lambda x0, y0, z0: 5*(z0-15)/r(x0, y0, z0)
	polymer.fix(1, "polymer nve")
	polymer.fix(2, "polymer langevin", t, t, 1.5, np.random.randint(2, high = 200000))
	polymer.fix(3, "polymer spring tether", 10, i/2, "NULL NULL", 0)
	polymer.timestep(0.02)
	polymer.compute("com polymer com")
	polymer.variable("ftotal equal fcm(polymer,x)")
	polymer.thermo_style("custom v_ftotal")
	polymer.thermo(1)
	print(i)
	polymer.run(5000)
	print(i)
	l = polymer.runs[0][0][0][1000:] + [i/2]
	u = [np.mean(l), i/2]
	l = polymer.extract_global("c_com", 1)
	print(i)
	np.savetxt("trial%dmean.txt" % i, u)
	np.savetxt("trial%dall.txt" % i, l)
	return u

data(15)
# num_cores = multiprocessing.cpu_count()

# results = Parallel(n_jobs = num_cores)(delayed(data)(i) for i in range(1, 90))
# Initialize random positions of particles
def position_init(length = 10):
	x1 = random.uniform(0, 28)
	y1 = random.uniform(0, 30)
	z1 = random.uniform(0, 30)
	i = 1
	Lx, Ly, Lz = [x1], [y1], [z1]
	while i < length:
		xi = random.uniform(Lx[i - 1] - 1, Lx[i - 1] + 1)
		if xi in Lx:
			while xi in Lx:
				xi = random.uniform(Lx[i - 1] - 1, Lx[i - 1] + 1)
		yi = random.uniform(Ly[i - 1] - 1, Ly[i - 1] + 1)
		if yi in Ly:
			while yi in Ly:
				yi = random.uniform(Ly[i - 1] - 1, Ly[i - 1] + 1)
		zi = random.uniform(Lz[i - 1] - 1, Lz[i - 1] + 1)
		if zi in Lz:
			while xi in Lx:
				xi = random.uniform(Lx[i - 1] - 1, Lx[i - 1] + 1)
		r = np.sqrt((Lx[i - 1]-xi)**2 + (Ly[i - 1]-yi)**2 + (Lz[i - 1]-zi)**2)
		xi = Lx[i - 1] + (Lx[i - 1] - xi)/r
		yi = Ly[i - 1] + (Ly[i - 1] - yi)/r
		zi = Lz[i - 1] + (Lz[i - 1] - zi)/r
	return zip(Lx, Ly, Lz)

# Eliminate repeated points and convert them to normal coordinate
def position_data_modify(x, y, z):
	i = 0
	j = 1
	n = 0
	while i < len(x):
	   while j < len(x):
	      if x[i] == x[j]:
	         if y[i] == y[j] and z[i] == z[j]:
	            x = x[:j]+x[j+1:]
	            y = y[:j]+y[j+1:]
	            z = z[:j]+z[j+1:]
	            j = j - 1
	      j += 1
	   i += 1
	   j = i+1

def CheckRepeats(x, y, z):
	i = 0
	j = 1
	n = 0
	while i < len(x):
	   while j < len(x):
	      if x[i] == x[j]:
	         if y[i] == y[j] and z[i] == z[j]:
	            return False
	      j += 1
	   i += 1
	   j = i+1
	return True