# Target gas mass vs. distance from black hole
# Showing how much the gas has refined
# Somewhat out of date; current refinement uses target gas radios instead
# Using DataLoader (python 3)

# Imports
from readData.DataLoader import DataLoader
#import h5py
#import simread.readsnapHDF5 as rs
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import math

print("Starting")
# Directory path and snapshot
directory = '/orange/paul.torrey/rlosacco/Spatial_Refine/'
version = 'vB'  # I used a weird naming scheme, don't worry about it
snapnum = 3
redshift = 0

# Simulation info
#snapname = "snapdir_"+snapnum+"/snapshot_"+snapnum
NumFilesPerSnapshot = 16
boxsize = 20000.
path = directory + version + '/'# + snapname
h = 0.7                        # Change h later
volume = (20./h)**3.
sim_mass_unit = 1e10/h
tgm_max = 1e-3 * sim_mass_unit
xaxis_max = 10

# Exponential refinement function
def expon(refine_factor, x, x_max, x_min):
  tgm_min = tgm_max / refine_factor
  coeff = np.power(refine_factor, (1/(x_min - x_max)))
  y = tgm_max * np.power(coeff, (x_max - x))
  return tgm_min, y

# v2b:
xmin = 1.0
xmax = 5.0
func_x = np.linspace(xmin, xmax, 50)
tgm_min, func_y = expon(10.0, func_x, xmax, xmin)

# Reading data
print("Reading data")
gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['BH_Distance', 'Masses'])
D2BH = gascat['BH_Distance']
gas_mass = gascat['Masses'] * sim_mass_unit
#D2BH = rs.read_block(path, "D2BH", parttype=0)    # For some reason this doesn't always work...
#gas_mass = rs.read_block(path, "MASS", parttype=0) * sim_mass_unit

# Use log10
loggasmass = np.log10(gas_mass)
tgm_min = math.log10(tgm_min)
tgm_max = math.log10(tgm_max)
func_y = np.log10(func_y)

# Plot
print("Plotting")
fig, ax = plt.subplots(figsize=(5,4))
ax.hist2d(D2BH, loggasmass, bins=128, range=[[0, xaxis_max], [5,8]])
# Vertical dashed lines
ax.axvline(x=xmin, ls='--', c='w', alpha=0.5)
ax.axvline(x=xmax, ls='--', c='w', alpha=0.5)
# Horizontal lines
print(tgm_min)
ax.axhline(y=tgm_min, xmin=0, xmax=xmin/xaxis_max, c='w', alpha=0.5) 
ax.axhline(y=tgm_max, xmin=xmax/xaxis_max, xmax=1, c='w', alpha=0.5)           
# function overplot
ax.plot(func_x, func_y, c='w', alpha=0.5)
ax.set_xlabel(r'$\mathrm{Distance\;[kpc]}$')
ax.set_ylabel(r'$\mathrm{log(Gas\;Cell\;Mass)\;}[\mathrm{M}_\odot]$')
fig.subplots_adjust(left=0.14, bottom=0.13, right=0.98, top=0.98)
fig.savefig('refine_'+version+'z'+str(redshift)+'.png')
plt.close()
