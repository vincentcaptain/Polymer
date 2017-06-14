from __future__ import division
from pylab import *
import os

# Time series of some correlation function plot and data saved. 1st argument is the raw data file, and
# the second one is the variable name you want to save as. 
o = sys.argv[1]
name = sys.argv[2]
da = loadtxt(o)
# Raw value of function f(0)f(t)
f0 = da[:,0][0]
ft = da[:,0]
f0t = []
f0t_ave = []
for i in range(len(da[:,0])):
	g0t = ft[i]*f0
	f0t.append(g0t)

for j in range(len(da[:,0])):
	if j == 0:
		f0t_ave.append(f0t[j])
	else:
		g0t_ave = (f0t_ave[j - 1] * j + f0t[j]) / (j + 1)
		f0t_ave.append(g0t_ave)
savetxt("correlation.txt", f0t_ave)
plot(range(len(f0t_ave)), f0t_ave)
xlabel("time")
ylabel("correlation_" + name)
savefig(name)
sys.exit()
