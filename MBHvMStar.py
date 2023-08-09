# Black hole mass vs. stellar mass of galaxies
# Note: can be combined with MBHvVelDisp.py for efficiency
# Recreating figures from Sijacki+15
# Fig 4: log(M_BH) vs. log(stellar halo mass) at z=0
# Fig 5: same at z=4,3,2,&1
# Using simread from torreylabtools (python 2.7)

# Imports
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,'/home/rlosacco/torreylabtools/Python')  # change this path to your directory
import simread.readsubfHDF5 as rsub
import simread.readhaloHDF5 as rhalo

print "Starting"
# Directory path and snapshot
version = 'debugging_level7'
path = '/orange/paul.torrey/rlosacco/'+version+'/'
snapnum = 171
a = 0.99881612
h = 0.6909
box = 20.

z = (1./a) - 1.
volume = (box/h)**3.
sim_mass_unit = 1e10/h

# Histogram Info
low_err = 10                                              # percentile for lower error
high_err = 90                                             # percentile for upper error
histrange = [(8., 13.), (5, 11)]                         # range in the y-direction
statsrange = (8., 13.)
histbins = 50                                             # number of bins
binsize = (statsrange[1] - statsrange[0])/histbins        # bin size

# Read in data
meds = []      # median
binmids = []   # middle of bin
pyerr = []     # percentile errors in the y-direction

# Example: cat =  rsub.subfind_catalog(path, snapshot, keysel=[])
catsub = rsub.subfind_catalog(path, snapnum, keysel=['SubhaloMassInRadType', 'GroupNsubs'])    # stellar mass from subfind
stellarHM = catsub.SubhaloMassInRadType[:, 4] * sim_mass_unit
cathalo = rhalo.HaloReader(path, snapnum, snapbase='snapshot')                                 # BH from halo reader
BHHM = []
for groupn in range(len(catsub.GroupNsubs)):
  for halon in range(catsub.GroupNsubs[groupn]):
    halo_idx = sum(catsub.GroupNsubs[:groupn]) + halon        # halo no. from stellarHM
    BHMA = cathalo.read('BHMA',  5, groupn, halon)
    if BHMA is None:
      BHHM.append(0.)
    else:
      BHHM.append(np.sum(BHMA) * sim_mass_unit)

# Make numpy array
BHHM = np.asarray(BHHM)
# Use log10
logstellar = np.log10(stellarHM)
logBHHM = np.log10(BHHM)

print("Getting stats")
vmeds = []    # for this simulation (v)ersion
vlow = []
vhigh = []
for b in range(histbins):
  # bins edges
  leftedge = statsrange[0] + binsize*b
  rightedge = statsrange[0] + binsize*(b+1)
  yvalues = []
  for x in range(len(logstellar)):
    if logstellar[x] >= leftedge and logstellar[x]<rightedge:
      # y values in the given bin
      if BHHM[x] == 0.:                                          # if there are no black holes, set 1e-5
        continue
      else:
        yvalues.append(logBHHM[x])
  median = np.median(yvalues)
  if len(yvalues) == 0:                                          # if bin is empty, errors are NaN
    vlow.append(median)
    vhigh.append(median)
  else:
    plow = np.percentile(yvalues, low_err)
    phigh = np.percentile(yvalues, high_err)  
    vlow.append(plow)
    vhigh.append(phigh)
  vmeds.append(median)
meds.append(vmeds)                                               # this version's medians added to list for all versions
 
vbinmids = []                                                    # middle of the bin
for hb in range(histbins):
  leftedge = statsrange[0] + binsize*hb
  rightedge = statsrange[0] + binsize*(hb+1)
  binmiddle = (leftedge + rightedge)/2.
  vbinmids.append(binmiddle)
binmids.append(vbinmids)
# Make numpy arrays
vlow = np.asarray(vlow)
vhigh = np.asarray(vhigh)
vmeds = np.asarray(vmeds)
vyerr = (vmeds-vlow, vhigh-vmeds)                                # length of the error
pyerr.append(vyerr)

# Plotting
print "Plotting"
fig, ax = plt.subplots(figsize=(5,5))
ax.hist2d(logstellar, logBHHM, bins=histbins, range=histrange)   # just for this (v)ersion
ax.errorbar(vbinmids, vmeds, yerr=vyerr, c='w', alpha=0.5)              
ax.set_xlabel(r'$\log(\mathrm{M}_\mathrm{stellar})[\mathrm{M}_\odot]$')
ax.set_ylabel(r'$\log(\mathrm{M}_\mathrm{BH})[\mathrm{M}_\odot]$')
fig.savefig('MBHvMStar_snap'+str(snapnum)+version+'.png')
plt.close()
 
# Make numpy arrays
meds = np.asarray(meds)
binmids = np.asarray(binmids)

# Plotting binned
print('Plotting binned')
fig, ax = plt.subplots(figsize=(5,5))
for i in range(len(meds)):
  ax.errorbar(binmids[i], meds[i], yerr=pyerr[i])
ax.set_xlabel(r'$\log(\mathrm{M}_\mathrm{stellar})[\mathrm{M}_\odot]$')
ax.set_ylabel(r'$\log(\mathrm{M}_\mathrm{BH})[\mathrm{M}_\odot]$')
ax.set_xlim(histrange[0])
ax.set_ylim(histrange[1])
fig.savefig('MBHvMStar_snap'+str(snapnum)+'binned.png')
plt.close()
