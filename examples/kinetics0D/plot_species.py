#!/usr/bin/python
##!/bender1/data/shared/software/apps/anaconda/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

print "Starting plot generation..."

# read cleaned data
myfile =  sys.argv[1]
print "Reading file: " + myfile

# The passed in scaling dir basename
fname =  sys.argv[2]
print "From basedir: " + fname

#get data
data = np.genfromtxt(myfile, delimiter=',', dtype=float, skip_header=1)

# extract header
headerdata = np.genfromtxt(myfile, delimiter=',', dtype=str, names=True)
header = np.array(headerdata.dtype.names)

print
print "Header data successfully extracted from logs!"
print header

# sort the data by the first (timestep) column
data = data[data[:,0].argsort()]
#print data

# convenience for basic plots
n_species = len(data[0,:])
print "plotting data for Temp plus " + str(n_species - 2) + " species"

timestep   = data[:,0]

# basic np vs total time
fig, ax1 = plt.subplots()

ax1.set_xlabel('Timestep')
ax1.plot(timestep , data[:,1], 'r' , label="Temp")
ax1.set_ylabel('Temperature', color='r')
ax1.tick_params('y', colors='r')

ax2 = ax1.twinx()
ax2.set_ylabel("Species Mass Fraction")
#ax2.tick_params('y')
for i in range(2, n_species):
    ax2.plot(timestep,data[:,i], label=str(header[i]))


box = ax2.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left' , bbox_to_anchor=(1.15, 0.5))

ax2.set_title(fname)

fig.savefig(fname + '.png')
plt.clf()


print "Plots created successfully from ", myfile
