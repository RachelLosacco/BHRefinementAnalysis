# Tully-Fisher reation
# Velocity dispersion vs. stellar mass
# Using DataLoader (python 3)

# Imports
from readData.DataLoader import DataLoader
import numpy as np
import matplotlib.pyplot as plt

print('Starting')
# Directory path and snapshot
version = "debugging_level7"                              
path = "/orange/paul.torrey/rlosacco/"+version+"/"
snapnum = 171
a = 0.99881612
box = 20.  # Mpc^3
h = 0.6909

# Info
z = (1./a) - 1.
volume = (box/h)**3.
sim_mass_unit = 1e10/h
fig, ax = plt.subplots(figsize=(5.5,4.5))
xlimits = [1e9, 1e12]
ylimits = [10., 1e3]
logxlimits = [9., 12.]

# Histogram Info
low_err = 10                                              # percentile for lower error
high_err = 90                                             # percentile for upper error
histrange = [xlimits, ylimits]
histbins = 50                                             # number of bins
binsize = (logxlimits[1] - logxlimits[0])/histbins              # bin size

# Data Info
parttypes = [4]
#keys = ["SubhaloStellarPhotometricsRad", "SubhaloStellarPhotometricsMassInRad"]
keys = ["SubhaloVelDisp", "SubhaloMassType"]    # circ vel from group data, then total stellar halo mass

# Read data
cat = DataLoader(path, part_types=parttypes, snap_num=snapnum, keys=keys)   
velocity = cat["SubhaloVelDisp"]  # km/s
mass = cat["SubhaloMassType"][:,4] * sim_mass_unit

# Get medians
print("Getting stats")
meds = []      # median values, len = histbins
lows = []      # lower errors, higher errors
highs = []
binmids = []   # middle of the bin
for b in range(histbins):
  # bins edges and middle
  leftedge = 10**(logxlimits[0] + binsize*b)
  rightedge = 10**(logxlimits[0] + binsize*(b+1))
  binmiddle = (leftedge + rightedge)/2.
  binmids.append(binmiddle)
  yvalues = []
  for x in range(len(mass)):
    if mass[x] >= leftedge and mass[x]<rightedge:
      # y values in the given bin
      yvalues.append(velocity[x])
  median = np.median(yvalues)
  meds.append(median)
  if len(yvalues) == 0:                     # if bin is empty, errors are NaN
    lows.append(median)
    highs.append(median)
  else:
    plow = np.percentile(yvalues, low_err)
    phigh = np.percentile(yvalues, high_err)
    lows.append(plow)
    highs.append(phigh)

# Make numpy arrays
lows = np.asarray(lows)
highs = np.asarray(highs)
meds = np.asarray(meds)
yerrs = (meds-lows, highs-meds)              # distance from median
yerrs = np.asarray(yerrs)

# Plotting
print("Plotting")
ax.scatter(mass, velocity, s=2)
ax.errorbar(binmids, meds, yerr=yerrs)
ax.set_ylabel("Velocity [km/s]")
ax.set_xlabel("Stellar Mass [$M_{\odot}$]")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(ylimits)
ax.set_xlim(xlimits)
fig.savefig('TullyFisher_snap'+str(snapnum)+version+'.png')
plt.close()