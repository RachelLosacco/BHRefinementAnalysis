# How long has my simulation been running? What are the lbt of each snapshot?
# Get scale factor from snapshot data
# Using DataLoader (python 3)

# Imports
import numpy as np
from scipy import interpolate
from units.springel_units import age_from_a
from readData.DataLoader import DataLoader
import h5py

# Path and snapshots
version = 'debugging_level7'
snapnums = np.arange(0, 171, 1)      # range of snapshots
path = '/orange/paul.torrey/rlosacco/'+version+'/'

# Get scale factors
print("Starting")
scalefactors = []
for i in snapnums:
    file_name = path+'snapdir_'+str(i).zfill(3)+'/snapshot_'+str(i).zfill(3)+'.0.hdf5'
    with h5py.File(file_name, "r") as ofile:
        header = ofile['Header']
        a = header.attrs['Time']
        scalefactors.append(a)

# Cosmo. parameters (from param.txt)
H0=0.6909
Om0=0.301712
Ode0=0.698288

# Calculate
scalefactors = np.asarray(scalefactors)
print('scale factors: ', scalefactors)
lbts = age_from_a(scalefactors, H0=H0, Om0=Om0, Ode0=Ode0)
print('lbt in gyr: ', lbts)
ages = age_from_a(0., H0=H0, Om0=Om0, Ode0=Ode0) - lbts
print('age of simulation in gyr: ', ages)

print('This simulation has been refining for ', max(ages)-min(ages), 'Gyr')