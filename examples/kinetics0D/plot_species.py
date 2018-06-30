#!/usr/bin/python

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

# The timestep size
dt =  float(sys.argv[3])
print "with dt: " + str(dt)

#get data
data = np.genfromtxt(myfile, delimiter=',', dtype=float, skip_header=1)

# extract header
headerdata = np.genfromtxt(myfile, delimiter=' ', dtype=str, names=True)
header = np.array(headerdata.dtype.names)
print
print "Header data successfully extracted from logs!"
print header

# no need to plot nonreacting species
species_dilutants = ["Ar", "AR", "H2", "O2", "H2O"]
#species_dilutants = ["Ar", "AR"]
exclude_species_idx = [0,1]; #also exclude time,Temperature data
for h in header:
    for sd in species_dilutants:
        if ( h == sd ):
            i, = np.where( header == h )
            exclude_species_idx.extend(i)
#print "Excluding from species data: " + str(header[exclude_species_idx])
all_indicies = [i for i in xrange(header.size)]
ind = [x for x in all_indicies if (x not in exclude_species_idx)]
print "plotting species: " + str(header[ind])

#get the time data
time = data[:,0]*dt

#get handles to figure and axis
fig, ax1 = plt.subplots()

# plot the Temperature. assume its always at column 1
ax1.set_xlabel('Time')
ax1.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#ax1.ticklabel_format(useOffset=False) # dont use +1e3 for temps
ax1.plot(time, data[:,1] , 'r' , label="Temp")
ax1.set_ylabel('Temperature', color='r')
ax1.tick_params('y', colors='r')

#double the y axis to plot temp alongside species with differnt scales
ax2 = ax1.twinx()
ax2.set_ylabel("Species Mass Fraction")

# Plot the species data
for i in ind:
    ax2.plot( time, data[:,i], label=str(header[i]) )

# put legend outside of plot center right
box = ax2.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(loc='center left' , bbox_to_anchor=(1.15, 0.5))

#title and save
ax2.set_title(fname.title() + " Combustion")
fig.savefig(fname + '.png')
plt.clf()

print "Plots created successfully from: ", myfile
