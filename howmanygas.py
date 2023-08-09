# Quickly see how many gas cells are near the black holes
# Creates histogram of gas particles vs. distance to nearest BH
# y-axis can be count, gas mass, gas radius, or density
# Un-comment lines in plotting section to choose
# Using DataLoader (python 3)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader

# Directory path and snapshot
version = 'Temporal_20_6'
path = '/orange/paul.torrey/rlosacco/'+version+'/'
snapnum = 796
galaxies = range(50)

# This is also in the AREPO code
def get_minradius(scalefactor):
  R_n = [300., 100., 50., 10., 5., 1., 0.5, 0.1]
  a_n = [0.72506334, 0.92045631, 0.99022114, 0.99735899, 0.99878656, 0.99950035, 0.99985724, 0.99992862, 1.0]
  for i in range(len(R_n)):
    if a_n[i] < scalefactor and scalefactor < a_n[i+1]:
      return R_n[i]
    else:
      print("could not determine minimum radius")

# Other info
fig, ax = plt.subplots(figsize=(5,4))
xaxis_max = 1.0
xaxis_min = -2.0
histrange = (-2., 3.)
histbins = 100
deltabin = (histrange[1] - histrange[0]) / (histbins + 1)   # size of each bin, in log, float
deltav = []                                           # an array of 50 with the volume of each shell, not in log
for shell in range(histbins):
  logadding_upr = histrange[0] + (shell+1)*deltabin   # from the start to the far edge in log
  adding_upr = 10**logadding_upr                      # out of log
  adding_upv = 4./3. * np.pi * (adding_upr**3.)       # volume of large sphere
  logsubtractingr =  histrange[0] + shell*deltabin    # from the start to the inner edge in log
  subtractingr = 10**logsubtractingr                  # out of log
  subtractingv = 4./3. * np.pi * (subtractingr**3.)   # volume of the smaller sphere
  deltav.append(adding_upv - subtractingv)            # the shell is the difference between the sphere volumes

# Reading in Data
for n in galaxies:
  print(n)
  #print(path)
  gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['Coordinates', 'Masses', 'Density'], sub_idx=0, fof_idx=n)
  bhcat = DataLoader(path, part_types=[5], snap_num=snapnum, keys=['Coordinates'], sub_idx=0, fof_idx=n)
  a = gascat.time
  #minradius = get_minradius(a)
  h = gascat.h
  sim_mass_unit = 1e10/h
  gaspos = gascat['Coordinates']
  gasmass = gascat['Masses'] * sim_mass_unit
  gasrho = gascat['Density'] * sim_mass_unit * h**3.
  bhpos = bhcat['Coordinates']
  
  # Not all snapshots & haloes have black holes
  if len(bhpos) == 0:
    print('no bh')
    continue
  else:    
    gasx = gaspos[:,0]
    gasy = gaspos[:,1]
    gasz = gaspos[:,2]
    bhx = bhpos[0,0]
    bhy = bhpos[0,1]
    bhz = bhpos[0,2]
    dx = gasx - bhx
    dy = gasy - bhy
    dz = gasz - bhz
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    mass_array, mass_binedges = np.histogram(np.log10(dist), bins=histbins, weights=gasmass, range=histrange)
    dist_array, dist_binedges = np.histogram(np.log10(dist), bins=histbins, range=histrange)
    
    gassize = ((3. * gasmass)/(4*np.pi* gasrho))**(1./3.)
    
    # Plotting histogram
    #plt.axhline(minradius, 0, 1, linestyle='--', color='k', alpha=0.5)    # size
    #ax.bar(dist_binedges[:-1], mass_array/deltav, width=deltabin, align='edge')                  # density
    #ax.hist(np.log10(dist), bins=histbins, range=histrange, weights=gasmass, histtype='step')    # mass
    #ax.hist(np.log10(dist), bins=histbins, range=histrange, histtype='step')                     # count
    ax.hist(np.log10(dist), bins=histbins, range=histrange, weights=gassize, histtype='step')    # size
ax.set_xlabel('log(Distance [kpc])')
#ax.set_ylabel('Gas Density [$M_\odot$/kpc$^3$]')            # density
#ax.set_ylabel('Gas Mass [$M_\odot$]')                       # mass
ax.set_ylabel('Gas Radius [kpc]')                           # size
ax.set_yscale('log')
ax.set_xlim(xaxis_min, xaxis_max)     # x-axis is already log
#ax.set_ylim(0, 8)           # size, y-linear
ax.axvline(1., 0, 1, color='grey', linestyle='--')
#fig.savefig('howmanyrho_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')           # density
#fig.savefig('howmanyMgas_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')          # mass
#fig.savefig('howmanygas_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')           # count
fig.savefig('howmanyRgas_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')          # size
plt.close()

