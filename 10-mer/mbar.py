#!/usr/bin/env python
from __future__ import division
from pylab import *
import os

output=open('metafile.txt','w')

betarange=arange(26, 34, 0.5)

for i in betarange:
	b = np.loadtxt("%d.txt" % i)
	y = [b[j][1] for j in range(int(len(b)/2), len(b))]
	x = list(range(len(y)))
	t_s = zip(x, y)
	print(len(x))
	savetxt('time_%d.txt' %i,t_s)
  	print >> output, 'time_%d.txt' %i, i, 1, len(x)

print "Run MBAR calculation"
print "syntax: ./mbar metafile weightsfile minQ maxQ bins biastype=(1,2)"

output.close()

os.system("./mbar metafile.txt weights.txt 23 37 25 2")
l1 = []
l2 = []
for x in betarange:
	b = np.loadtxt("%s.txt" % x)
	k = [b[i][1] for i in range(int(len(b)/2), len(b))]
	h = np.histogram(k, bins = 100)
	# p = [-np.log(h[0][i]/np.sum(h[0])) for i in range(len(h[0]))]
	pu = -(np.log(h[0][:]/np.sum(h[0]))+0.5*(h[1][1:]-0.5*(h[1][1]-h[1][0])-x)**2)
	l, = plt.plot(h[1][1:], pu)
	l1.append(l)
	l2.append("%s" % x)
	print(x)



plt.legend(l1, l2)
plt.xlabel("x")
plt.ylabel("unbiased potential")
plt.savefig("unbiased")
plt.show()

data = loadtxt("phi.txt")
plt.plot([data[i][0] for i in range(len(data))], [data[i][1] for i in range(len(data))])
plt.xlabel("x")
plt.ylabel("Pu")
plt.savefig("mbar")
plt.show()
