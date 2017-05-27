#!/usr/bin/env python
from __future__ import division
from pylab import *
import os

output=open('metafile.txt','w')

betarange=arange(23,38,1)

for i in betarange:
	b = np.loadtxt("time_%s.txt" % i)
	y = [b[j][1] for j in range(len(b)) if j % 10 == 0]
	x = [b[j][0] for j in range(len(b)) if j % 10 == 0]
	t_s = zip(x, y)
	print(len(x))
	np.savetxt("time_%s.txt" % i, t_s, fmt = "%f")
	#print >> output, 'time_%s.txt' %i, 1, len(b)
	output.write("time_%s.txt %s 1 %d" %(i, i, len(y)))
	print "Creating input files"

print "Run MBAR calculation"
print "syntax: ./mbar metafile weightsfile minQ maxQ bins biastype=(1,2)"

output.close()

os.system("./mbar metafile.txt weights.txt 15 45 50 2")
l1 = []
l2 = []
for x in betarange:
	b = np.loadtxt("%s.txt" % x)
	k = [b[i][1] for i in range(int(len(b)/2), len(b))]
	h = np.histogram(k, bins = 100)
	p = [h[0][i]/np.sum(h[0]) for i in range(len(h[0]))]
	#pu = [np.log(i)+5*(h[-x)**2 for i in p]
	pu = np.log(h[0][:]/np.sum(h[0]))+0.5*(h[1][1:]-0.5*(h[1][1]-h[1][0])-x)**2
	l, = plt.plot(h[1][1:], pu)
	l1.append(l)
	l2.append("%s" % x)
	print(x)



plt.legend(l1, l2)
plt.xlabel("x")
plt.ylabel("unbiased")
plt.savefig("unbiased by hand")
plt.show()