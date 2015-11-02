#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import glob, os

os.chdir("plotfiles/")

bdvplots=[]
for file in glob.glob("bdvline*"):
    bdvplots.append(file)

asplots=[]
for file in glob.glob("asline*"):
    asplots.append(file)

blplots=[]
for file in glob.glob("blline*"):
    blplots.append(file)

asplots.sort()
blplots.sort()
bdvplots.sort()
	
print bdvplots
print asplots
print blplots

for i in range(0,len(blplots)):
	with open(blplots[i]) as f:
	    data = f.read()
	data = data.split('\n')[:-1]
	
	x = [float(row.split()[0]) for row in data]
	y = [float(row.split()[1]) for row in data]
	
	plt.title("Plot_"+str(i))    
	plt.xlabel('share')
	plt.ylabel('stretch')

	plt.plot(x,y,'--s', ms=10, markerfacecolor="None", markeredgecolor='blue', markeredgewidth=2)
	
	with open(asplots[i]) as f:
	    data = f.read()
	data = data.split('\n')[:-1]
	
	x = [float(row.split()[0]) for row in data]
	y = [float(row.split()[1]) for row in data]
		
	plt.plot(x,y,'^', ms=10, markerfacecolor="None", markeredgecolor='green', markeredgewidth=2)
	
	with open(bdvplots[i]) as f:
	    data = f.read()
	data = data.split('\n')[:-1]
	
	x = [float(row.split()[0]) for row in data]
	y = [float(row.split()[1]) for row in data]
		
	plt.plot(x,y,'o', ms=10, markerfacecolor="None", markeredgecolor='red', markeredgewidth=2)

	plt.savefig('plot'+str(i)+'.png', bbox_inches='tight')
	
	plt.close()
	
