#!/usr/bin/env python

from pylab import *
import os

print "Creating gaussian dataset"

output=open('metafile.txt','w')
numpoints=5000
betarange=arange(-5,6,1)

for i in betarange:
  bb=np.random.randn(numpoints)+i
  aa=arange(numpoints)
  zz=zip(aa,bb)
  savetxt('time_%d.txt' %i,zz)
  print >> output, 'time_%d.txt' %i, i, numpoints

print "Creating input files"

print "Run MBAR calculation"
print "syntax: ./mbar metafile weightsfile minQ maxQ bins biastype=(1,2)"

output.close()

os.system("./mbar metafile.txt weights.txt -8 8 80 1")

da=loadtxt('phi.txt')
da1=loadtxt('psi.txt')

subplot(121)
plot(da[:,0],da[:,1],'ob')
x=arange(-8,8,.1)
plot(x,-x**2./2.,'-r')
xlabel(r'$Q$')
ylabel(r'$\phi(Q)$')

subplot(122)
plot(da1[:,0],da1[:,1],'ob')
plot(betarange,betarange**2./2.,'xr',lw=4,mew=2)
xlabel(r'$\beta$')
ylabel(r'$\psi(\beta)$')

savefig('test.pdf')
