from joblib import Parallel, delayed
import multiprocessing
from pylab import *
import os

omega = arange(-5, 5, 0.1)

def run_FFS(o):
	os.system("python2 FFS.py %s" % o)

Parallel(n_jobs=100)(delayed(run_FFS)(i) for i in omega)
