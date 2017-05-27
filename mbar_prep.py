import os
import string
from math import *
from numpy import *
import sys

home=os.getcwd()

# Gets list of umbrella directories
folder=os.listdir('%s/' %home)

k=1
T=eval(sys.argv[1])
equil=6000
subsample=50
output=open("metafile.txt","w")

for i in folder:
	try: 
		temp=os.path.getsize('%s/sample.txt' %i)
	except:
		continue
		
	if temp==0: os.system('rm %s/sample.txt' %i)
	
	try:
		da=loadtxt('%s/sample.txt' %i)
	except:
		continue

	if len(da[:,0])!=0:
		
		print i
		#os.system("grep 'variable' %s/prod.in > %s/umbrella.txt" %(i,i))
		# Time
		t=da[equil:,0]
		# Biased Order Parameter
		q=da[equil:,5]
		# Binned Order Parameter
		q2=da[equil:,6]

		step=len(q)/subsample

		# Creates a restricted time_series of just pertinent data
		# Useful construct for sub sampling
		
		output_time=open('%s/time_series.txt' %i,'w')
		for j in range(step): print >> output_time, t[subsample*j], q[subsample*j], q2[subsample*j]
		output_time.close()
		
			
		input= open("%s/umbrella.txt" %i).read()
		inp=[string.split(line," ") for line in string.split(input,"\n")]
		
		# Umbrella Potential Minima
		q_min=eval(inp[0][3])
		
		# Spring Constant in Units of kT, MBAR program assumes k=1
		
		print '%f' %q_min, '%f' %(kappa), step
		# Metadata file
		print >> output, i+'/time_series.txt', '%f' %q_min, '%f' %(kappa), step

	else: print '%s has a blank sample.txt' %i
 


# umbrella.txt has a format like:
# variable x index 39 <- order parameter minima

# sample.txt has a format like:
# time temperature pote_energy volume total_energy other biased_op desired_op
