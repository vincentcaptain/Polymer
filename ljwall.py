import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpi4py import MPI
me = MPI.COMM_WORLD.Get_size()
from lammps import lammps
lmp = lammps()

