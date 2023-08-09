# Comparing MBH of refined and non-refine (v0) runs
# Using simread from torreylabtools (python 2.7)

# Imports
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,'/home/rlosacco/torreylabtools/Python')
import simread.readsubfHDF5 as rsub
import simread.readhaloHDF5 as rhalo
import simread.readsnapHDF5 as rsnap

print "Starting"
cumulative = True  # cumulative histogram
large = True       # Includes larger range
# Simulation
version = 'debugging_level7'
path = '/orange/paul.torrey/rlosacco/'+version+'/'
snapnum = 171
a = 0.99881612
h = 0.6909
box = 20.

# Info
z = (1./a) - 1.
volume = (box/h)**3.
sim_mass_unit = 1e10/h
fig, ax = plt.subplots(figsize=(5,5))
if cumulative == True:
  cumout = '_cumulative'
  cumulative = -1
else:
  cumout = ''

# Histogram Info
low_err = 10                                              # percentile for lower error
high_err = 90                                             # percentile for upper error
if large == True:
  histrange = (7., 10.)
  largeout = 'large'
else:
  histrange = (6.0, 10.)
  largeout = ''
histbins = 50                                             # number of bins
binsize = (histrange[1] - histrange[0])/histbins        # bin size
meds = []      # median
binmids = []   # middle of bin
pyerr = []     # percentile errors in the y-direction

# Data Info
catsub = rsub.subfind_catalog(path, snapnum, keysel=['GroupNsubs', 'GroupLenType', 'SubhaloLenType'])

BHHM = []  # black hole halo mass
doubleBHcount = 0
# load in BHMA from all groups & subs
BHMA = rsnap.read_block(path+"snapdir_"+str(snapnum).zfill(3)+"/snapshot_"+str(snapnum).zfill(3), "BHMA", parttype=5)  
BH_idx = 0         # counting all BH
BH_groupidx = 0    # BH in group but not in subhalo
subhalo_idx = 0    # counting subhalos
for groupn in range(len(catsub.GroupNsubs)):
  BH_idx = BH_groupidx                                      # include BH not in subhalo in BH count
  for subhalon in range(catsub.GroupNsubs[groupn]):
    nBHs = catsub.SubhaloLenType[subhalo_idx+subhalon, 5]   # no. of BH in this subhalo
    if nBHs > 0:
      BH_mass = BHMA[BH_idx:nBHs+BH_idx]                    # BH mass array for BH in this subhalo
      BH_idx += nBHs
      BHHM.append(np.sum(BH_mass) * sim_mass_unit)          # sum all BH mass in subhalo
    else:
      BHHM.append(0.)
  BH_groupidx += catsub.GroupLenType[groupn, 5]             # all BH in group, to idx of next group
  subhalo_idx += catsub.GroupNsubs[groupn]                  # all subs in group, to idx of next group

# Make numpy array
BHHM = np.asarray(BHHM)
# Use log10
logBHHM = np.log10(BHHM)

'''
# Commented out because it didn't always work
# Error bars
vmeds = []
vlow = []
vhigh = []
for b in range(histbins):
  # bin edges
  leftedge = histrange[0] + binsize*b
  rightedge = histrange[0] + binsize*(b+1)
  yvalues = []
  for x in range(len(logBHHM)):
    if logBHHM[x] >= leftedge and logBHHM[x] <=rightedge:
      if BHHM[x] == 0.:
        continue
      else:
        yvalues.append(logBHHM[x])
  median = np.median(yvalues)
  if len(yvalues) == 0:
    vlow.append(median)
    vhigh.append(median)
  else:
    plow = np.percentile(yvalues, low_err)
    phigh = np.percentile(yvalues, high_err)
    vlow.append(plow)
    vhigh.append(phigh)
  vmeds.append(median)
meds.append(vmeds)

vbinmids = []
for hb in range(histbins):
  leftedge = histrange[0] + binsize*hb
  rightedge = histrange[0] + binsize*(hb+1)
  binmiddle = (leftedge + rightedge) / 2.
  vbinmids.append(binmiddle)
binmids.append(vbinmids)
vlow = np.asarray(vlow)
vhigh = np.asarray(vhigh)
vmeds = np.asarray(vmeds)
vyerr = (vmeds-vlow, vhigh-vmeds)
pyerr.append(vyerr)
'''

# Plotting
print("Plotting")
ax.hist(logBHHM, bins=histbins, range=histrange, cumulative=cumulative, histtype='step')
ax.set_xlabel(r'$\log(\mathrm{M}_\mathrm{BH})[\mathrm{M}_\odot]$')
fig.savefig('MBHcomp'+cumout+largeout+str(snapnum)+version+'.png')
plt.close()

