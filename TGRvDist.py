# Target gas radius vs distance from the nearest black hole
# Like TGMvDist.py and SizevDist.py
# Using DataLoader (python 3)

# Imports
from readData.DataLoader import DataLoader
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import math
from astropy import units as u
from astropy import constants as con

print("Starting")
# Directory path and snapshot
path = '/orange/paul.torrey/rlosacco/Temporal_20_7/'
snapnum = 173
galaxies = range(0, 200)
min_size = 0.3               
outname = 'TGRvDist_Temporal20x7_snap'+str(snapnum)+'.png'

# Simulation info
fig, ax = plt.subplots(figsize=(5,4))
xaxis_max = 1.5
yaxis_max = 1.5

# Reading data
for n in galaxies:
  print(n)  
  gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['Coordinates', 'Masses', 'Density'], sub_idx=0, fof_idx=n)
  bhcat = DataLoader(path, part_types=[5], snap_num=snapnum, keys=['Coordinates'], sub_idx=0, fof_idx=n)
  h = gascat.h
  a = gascat.time
  sim_mass_unit = 1e10*u.Msun / h
  sim_pos_unit = u.kpc / h
  sim_rho_unit = sim_mass_unit / sim_pos_unit**3.
  gaspos = gascat['Coordinates'] * sim_pos_unit
  gasmass = gascat['Masses'] * sim_mass_unit
  gasrho = gascat['Density'] * sim_rho_unit
  bhpos = bhcat['Coordinates'] * sim_pos_unit
  #bhdist = gascat['BH_Distance'] * sim_pos_unit   # Doesn't always work for some reason...
    
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
    dist = dist.decompose().to('kpc').value
  
  # Cell radius from mass and density
  gas_r = ((3.*gasmass)/(4.*np.pi*gasrho))**(1./3.)
  gas_r = gas_r.decompose().to('kpc').value
  
  ax.hist2d(dist, gas_r, bins=128, range=[[0, xaxis_max], [0,yaxis_max]])

# Refinement function
refine_factor = 3.0
x = np.linspace(0, 10, 100)   # from 0 to 10 kpc out from center
y = x/refine_factor

# Plot
print("Plotting")
ax.plot(x, y, color='white')
ax.axhline(min_size, 0, 1, color='yellow', linestyle='--', alpha=0.5)
ax.axvline(1., 0, 1, color='yellow', linestyle='--', alpha=0.5)
ax.set_xlabel(r'$\mathrm{Distance\;[kpc]}$')
ax.set_ylabel(r'$\mathrm{Gas\;Cell\;Radius\;[kpc]}$')
fig.subplots_adjust(left=0.14, bottom=0.13, right=0.98, top=0.98)
fig.savefig(outname)
plt.close()
