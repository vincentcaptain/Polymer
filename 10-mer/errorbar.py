from pylab import *
ax = subplot(111)

x=arange(-2,2,0.01)
f=-x**2.+x**4.
err=zeros(len(x))+0.5

plot(x,f,'-b',lw=4)

ft=f+3*err
fb=f-3*err

ftot=append(ft,fb[::-1])
xtot=append(x,x[::-1])

fill(xtot,ftot, facecolor='0.5', edgecolor='None',alpha=0.25)

ft=f+2*err
fb=f-2*err

ftot=append(ft,fb[::-1])
xtot=append(x,x[::-1])

fill(xtot,ftot, facecolor='0.5', edgecolor='None',alpha=0.25)

ft=f+err
fb=f-err

ftot=append(ft,fb[::-1])
xtot=append(x,x[::-1])

fill(xtot,ftot, facecolor='0.5', edgecolor='None',alpha=0.25)


xlabel("t_obs")
ylabel("msd")

show()
