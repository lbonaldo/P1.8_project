#!/usr/bin/env python

#to call c-routines
from ctypes import *
#include file
from math import *
#to get argument from terminal
from sys import argv

input_file = '../examples/argon_108.inp'

#import c-libraries
dso = CDLL("../sysdlib.so")

def output(sys_data, file1, file2):

    print("{} {:20.8f} {:20.8f} {:20.8f} {:20.8f}".format(sys_data.nfi, sys_data.temp, sys_data.ekin, sys_data.epot, sys_data.ekin+sys_data.epot))
    print("{} {:20.8f} {:20.8f} {:20.8f} {:20.8f}".format(sys_data.nfi, sys_data.temp, sys_data.ekin, sys_data.epot, sys_data.ekin+sys_data.epot), file=file1)
    print("{}\n nfi={} etot={:20.8f}\n".format(sys_data.natoms, sys_data.nfi, sys_data.ekin+sys_data.epot), file=file2)

    for i in range(sys_data.natoms):
        print("Ar  {:20.8f} {:20.8f} {:20.8f}\n".format(sys_data.rx[i], sys_data.ry[i], sys_data.rz[i]), file=file2)

def read_file(restfile, sys_data):
    rest = open(restfile, 'r')
    for i, line in enumerate(rest, 0):
        if (i < sys_data.natoms):
            sys_data.rx[i] = float(line.split()[0])
            sys_data.ry[i] = float(line.split()[1])
            sys_data.rz[i] = float(line.split()[2])
        elif (i >= sys_data.natoms and i < 2*sys_data.natoms):
            sys_data.vx[i-sys_data.natoms] = float(line.split()[0])
            sys_data.vy[i-sys_data.natoms] = float(line.split()[1])
            sys_data.vz[i-sys_data.natoms] = float(line.split()[2])
        else:
            print("Initilization error.")
    rest.close()

#structure with the constants
class mdsys(Structure):
    _fields_ = [("natoms", c_int), ("nfi", c_int), ("nsteps", c_int),
                ("dt", c_double), ("mass", c_double), ("epsilon", c_double),
                ("sigma",c_double), ("box",c_double), ("rcut",c_double),
                ("ekin",c_double), ("epot",c_double), ("temp",c_double),
                ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),
                ("vx", POINTER(c_double)), ("vy", POINTER(c_double)), ("vz", POINTER(c_double)),
                ("fx", POINTER(c_double)), ("fy", POINTER(c_double)), ("fz", POINTER(c_double))]

with open(input_file, 'r') as f:
    inputs = [i.split()[0] for i in f] ##add an error if not initialized? if(fp)?

#initializing the structure from file
sys_data = mdsys()

sys_data.natoms = int(inputs[0])
sys_data.mass = float(inputs[1])
sys_data.epsilon=float(inputs[2])
sys_data.sigma=float(inputs[3])
sys_data.rcut=float(inputs[4])
sys_data.box=float(inputs[5])
restfile=inputs[6]
trajfile=inputs[7]
ergfile=inputs[8]
sys_data.nsteps=int(inputs[9])
sys_data.dt=float(inputs[10])
nprint=int(inputs[11])

#"allocate" memory
sys_data.rx = (c_double * sys_data.natoms)()
sys_data.ry = (c_double * sys_data.natoms)()
sys_data.rz = (c_double * sys_data.natoms)()
sys_data.vx = (c_double * sys_data.natoms)()
sys_data.vy = (c_double * sys_data.natoms)()
sys_data.vz = (c_double * sys_data.natoms)()
sys_data.fx = (c_double * sys_data.natoms)()
sys_data.fy = (c_double * sys_data.natoms)()
sys_data.fz = (c_double * sys_data.natoms)()

#initialization
read_file(restfile, sys_data)
            
dso.azzero.argtypes = [ POINTER(c_double), c_int]
dso.azzero(sys_data.fx, sys_data.natoms)
dso.azzero(sys_data.fy, sys_data.natoms)
dso.azzero(sys_data.fz, sys_data.natoms)

sys_data.nfi = 0

# initialize forces and energies
dso.force(byref(sys_data))
dso.ekin(byref(sys_data))

erg = open(ergfile,"w")
traj = open(trajfile,"w")

print(f"Starting simulation with {sys_data.natoms} atoms for {sys_data.nsteps} steps.")
print("     NFI            TEMP            EKIN                 EPOT              ETOT")

output(sys_data, erg, traj)

# main MD loop
    
for sys_data.nfi in range(1, sys_data.nsteps + 1, 1):

    # write output, if requested
    if ((sys_data.nfi % nprint) == 0):
        output(sys_data, erg, traj)
        
    #  propagate system and recompute energies
    dso.velverlet1(byref(sys_data))
    dso.force(byref(sys_data))
    dso.velverlet2(byref(sys_data))
    dso.ekin(byref(sys_data))

# clean up: close files, free memory
print("Simulation Done.")

erg.close()
traj.close()
