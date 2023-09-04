#!/usr/bin/python

import sys
import math
import numpy as np

# read command line arguments, if any

# default values
mass_weighted = False
is_selection  = False

Nargs = len(sys.argv)

for n in range(1,Nargs):
    
    if sys.argv[n] in ("-h", "--help"):
        print ""
        print "usage: ", sys.argv[0], " [ --select n1 n2 ... ]"
        print ""
        print " options:"
        print "    --help          :    prints this help"
        print "    --mass-weighted :    calculates mass-weighted VACF (default = no)"
        print "    --select        :    select only atoms with indices n1 n2 ... (must be the last option) (default = all)"
        print ""
        exit()
        
    elif sys.argv[n] in ("-s", "--select"):
        is_selection = True
        Nselect = Nargs - n - 1
        selection = []
        selection = [ int(n)-1 for n in sys.argv[n+1:] ]    # -1 to make indices 0-based
    
    elif sys.argv[n] in ("-mw", "--mass-weighted"):
        mass_weighted = True


f = open('velocities.dat','r')

# first, get the number of frames and atoms
Nframes = 0
Natoms  = 0
atom_type = []

for line in f:

    if "frame" in line:
        Nframes += 1
    
    elif Nframes == 1:
        Natoms += 1
        column = line.split()
        atom_type.append(column[1])


# assign masses
mass = np.ones(Natoms)

if mass_weighted:
    for atom in range(Natoms):
        
        if atom_type[atom] in ("HA","HC","H"):
            mass[atom] = np.float64(1.008)
            
        elif atom_type[atom] in ("CA","CM","C"):
            mass[atom] = np.float64(12.011)
            
        elif atom_type[atom] in ("NA","N"):
            mass[atom] = np.float64(14.007)

        elif atom_type[atom] in ("Ox","O"):
            mass[atom] = np.float64(16.000)

        elif atom_type[atom] in ("Ti","TI"):
            mass[atom] = np.float64(47.867)



if not is_selection:
    Nselect = Natoms
    selection = range(Natoms)

# rewind file
f.seek(0)

print "# frames,atoms=", Nframes, Natoms, "   selected=",selection

v = np.zeros([Nframes,Natoms,3])
t = np.zeros(Nframes)
frame = -1

for line in f:
    
    if "#" in line:
        frame += 1
        column = line.split()
        t[frame] = np.float64(column[3])
    else:
        column = line.split()
        atom = int(column[0]) - 1      # index starts at zero, so the -1
        
        # mass == 1 if mass_weighted = False
        m = mass[atom]
        v[frame,atom,0] = np.float64(column[2])
        v[frame,atom,1] = np.float64(column[3])
        v[frame,atom,2] = np.float64(column[4])


# Calculate auto-correlation
ac = np.zeros(Nframes)

Nframe0 = 1
frame_jump = int(Nframes/Nframe0)

for frame0 in range(0,Nframes,frame_jump):
    
    ac0 = np.float64(0.0)
    v0 = v[frame0]
    
    for atom in selection:
        ac0 += mass[atom] * np.dot( v0[atom], v0[atom] )
    
    for frame in range(Nframes-frame0):
        
        for atom in selection:
            ac[frame+frame0] += mass[atom] * np.dot( v0[atom], v[frame,atom] ) / ac0


for frame in range(Nframes):
    print t[frame], ac[frame]/np.float64(Nframe0)       # wrong: normalization is different for each origin


