# Plotting stellar mass v halo mass
# Using simread from torreylabtools (python 2.7)

# Imports
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(1,'/home/rlosacco/torreylabtools/Python')

import simread.readsubfHDF5 as rsub

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
histrange = [(9., 14.), (7., 12.)]                         # range in the y-direction
statsrange = (9., 14.)
histbins = 50                                             # number of bins
binsize = (statsrange[1] - statsrange[0])/histbins        # bin size

# Observed model
def smhm_model(mass, afac):
  m0 = 11.59
  m1 = 1.195
  n0 = 0.0351
  n1 = -0.0247
  beta0 = 1.376
  beta1 = -0.826
  gamma0 = 0.608
  gamma1 = 0.329

  m = 10.**(m0 + (1 - afac) * m1) * 1e-10
  n = n0 + (1 - afac) * n1
  beta = beta0 + (1 - afac) * beta1
  gamma = gamma0 + (1 - afac) * gamma1

  val = 2. * n / ( (mass / m)**(-beta) + (mass / m)**gamma )

  return val * mass

# Binned medians
meds = []
binmids = []   # middle of bin
pyerr = []     # percentile errors in the y-direction

# Read in data
cat = rsub.subfind_catalog(path, snapnum, keysel=['SubhaloMassInRadType', 'SubhaloMass'])
Mstar = cat.SubhaloMassInRadType[:, 4] * sim_mass_unit
Mhalo = cat.SubhaloMass[:] * sim_mass_unit
logMstar = np.log10(Mstar)
logMhalo = np.log10(Mhalo)

# Make model data
mxrange = np.linspace(10**statsrange[0], 10**statsrange[1], 200)   # range of halo mass
model = smhm_model(mxrange/sim_mass_unit, a) * sim_mass_unit
logx = np.log10(mxrange)
logmodel = np.log10(model)

print("Getting stats")
vmeds = []    # for this simulation (v)ersion
vlow = []
vhigh = []
for b in range(histbins):
  # bins edges
  leftedge = statsrange[0] + binsize*b
  rightedge = statsrange[0] + binsize*(b+1)
  yvalues = []
  for x in range(len(logMhalo)):
    if logMhalo[x] >= leftedge and logMhalo[x]<rightedge:
      # y values in the given bin
      if Mstar[x] == 0.:                                          # if there are no black holes, set 1e-5
        continue
      else:
        yvalues.append(logMstar[x])
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

print('Plotting')
# Plotting
fig, ax = plt.subplots(figsize=(5,5))
ax.hist2d(logMhalo, logMstar, bins=histbins, range=histrange)
ax.errorbar(vbinmids, vmeds, yerr=vyerr, c='w', alpha=0.5, label="Median")              # just for this (v)ersion
ax.plot(logx, logmodel, 'y--', alpha=0.8, label="Moster+13")
ax.legend()
ax.set_xlabel(r'$\log(\mathrm{M}_\mathrm{halo})[\mathrm{M}_\odot]$')
ax.set_ylabel(r'$\log(\mathrm{M}_\mathrm{stellar})[\mathrm{M}_\odot]$')
fig.savefig('MstarvMhalo_snap'+str(snapnum)+version+'.png')
plt.close()

# Make numpy arrays
meds = np.asarray(meds)
binmids = np.asarray(binmids)

# Plotting binned
print('Plotting binned')
fig, ax = plt.subplots(figsize=(5,5))
for i in range(len(meds)):
  ax.errorbar(binmids[i], meds[i], yerr=pyerr[i])
ax.plot(logx, logmodel, 'k--', label="Moster+13")
ax.set_xlabel(r'$\log(\mathrm{M}_\mathrm{halo})[\mathrm{M}_\odot]$')
ax.set_ylabel(r'$\log(\mathrm{M}_\mathrm{stellar})[\mathrm{M}_\odot]$')
ax.set_ylim(histrange[1])
ax.legend()
fig.savefig('MstarvMhalo_snap'+str(snapnum)+'binned.png')
plt.close()