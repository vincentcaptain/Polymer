import numpy as np
import matplotlib.pyplot as plt

t = []
h = []
l = np.loadtxt("data.txt")
i = 0
ele = []
leg = []
j = 27
x = [l[j][1] for j in range(len(l))]
while i < len(x)/50001:
	t.append(x[i*50001:(i+1)*50001])
	i += 1

for x in t:
	h.append(np.histogram(x, bins = 100, normed = True))

for x in h[1:]:
	a, = plt.plot(x[1][1:], x[0])
	ele.append(a)
	leg.append("%i" % j)
	j += 0.5



plt.legend(ele, leg)
plt.show()