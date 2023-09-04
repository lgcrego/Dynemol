#!/usr/bin/python

import sys
import numpy as np


if( len(sys.argv) != 2 ):

    print ""
    print "Usage: ",sys.argv[0]," <.pdb file>"
    print ""
    sys.exit(1)

f = open(sys.argv[1],'r')


# first, get the number of atoms and species
Natoms = 0
species = []

for line in f:

    if "ATOM" in line:
        Natoms += 1
        species.append( line.split()[2] )
        
    elif "END" in line:
        break



# rewind file
f.seek(0)

frame = 0
t0 = np.float64(0)
r0 = np.zeros([Natoms,3])
r1 = np.zeros([Natoms,3])
v  = np.zeros([Natoms,3])


# now, calc. velocities
for line in f:
    
    if "TITLE" in line:
        frame += 1
        column = line.split()
        t1 = np.float64(column[5])
        
    elif "ATOM" in line:
        column = line.split()
        atom = int(column[1]) - 1      # index starts at zero, so the -1
        r1[atom] = [ np.float64(column[5]), np.float64(column[6]), np.float64(column[7]) ]
    
    elif "END" in line:
        if frame == 1:
            r0 = np.copy(r1)
            t0 = t1
            
        else:
            dt = t1 - t0
            v = (r1 - r0)/dt
            
            print "# frame,t,dt= ", frame, t1, dt
            
            for n in range(Natoms):
                print "  ",n+1, species[n], v[n][0], v[n][1], v[n][2]
            
            t0 = t1
            r0 = np.copy(r1)

f.close()
