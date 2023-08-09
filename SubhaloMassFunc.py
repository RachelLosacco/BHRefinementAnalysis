# Make subhalo stellar mass function
# Number of galaxies at a given stellar mass vs. stellar mass
# Using DataLoader (python 3)

# Imports
from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt

print('Starting')
# Directory path and snapshot
version = 'debugging_level7'
path = '/orange/paul.torrey/rlosacco/'+version+'/'
snapnum = 171
a = 0.99881612
box = 20.  # Mpc^3
h = 0.6909

# Info
z = (1./a) - 1.
volume = (box/h)**3.
sim_mass_unit = 1e10/h
fig, ax = plt.subplots(figsize=(5,5))

# Data Info
parttypes = [4]
keys = ["SubhaloMassType", "SubhaloGrNr"]
cat = DataLoader(path, part_types=parttypes, snap_num=snapnum, keys=keys)
subhalos = np.sort(cat["SubhaloMassType"][:,4] * sim_mass_unit)

# Plotting
print('Plotting')
ax.plot(subhalos[:-1][::-1], np.arange(len(subhalos)-1))
ax.set_ylabel("Cumulative Count")
ax.set_xlabel("Stellar Mass ($M_{\odot}$)")
ax.set_xscale("log")
ax.set_yscale("log")
fig.savefig('SubhaloMassFunc_snap'+str(snapnum)+version+'.png')
plt.close()
