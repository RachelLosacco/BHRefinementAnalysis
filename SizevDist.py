# Cell Size vs Distance to BH
# 2D histogram, very similar to TGRvDist.py
# Using DataLoader (python 3)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader


print('Starting')
# Change to your path
directory = '/orange/paul.torrey/rlosacco/debugging_level7'
output_name = './SizevDist_Plots/SizevDist_snap'+str(snapnum)+'.png' # where to save figure
box = 20.0

# Number of snapshots
NumFilesPerSnapshot = 16
path = directory # + snapname

# Reading data over range of snapshots
for snapnum in range(2, 171):
  print("Reading data in snap", snapnum)
  gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['Coordinates', 'Masses', 'Density'])
  bhcat = DataLoader(path, part_types=[5], snap_num=snapnum, keys=['Coordinates'])
  a = gascat.time
  h = gascat.h
  volume = (box/h)**3.
  sim_mass_unit = 1e10/h
  #n = 2  # fof index (central galaxy)
  #gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['BH_Distance', 'Masses', 'Density'], sub_idx=0, fof_idx=n) # single galaxy
  
  gaspos = gascat['Coordinates']
  bhpos = bhcat['Coordinates']
  # Not all haloes or snapshots have black holes
  if len(bhpos) == 0:
    print('no bh')
    D2BH = 1e10
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
    D2BH = np.sqrt(dx*dx + dy*dy + dz*dz)
  
  gas_mass = gascat['Masses'] * sim_mass_unit
  gas_rho = gascat['Density'] * sim_mass_unit
  gas_r = ((3.*gas_mass)/(4.*np.pi*gas_rho))**(1./3.)  # Cell radius from mass and density
  
  # Refinement function
  refine_factor = 3.0
  min_size = 3e-1               # minimum cell size 300pc = 0.3kpc
  x = np.linspace(0, 10, 100)   # from 0 to 10 kpc out from center
  y = x/refine_factor
  
  # Plotting
  print('Plotting')
  xaxis_max = 3
  yaxis_max = 3
  fig, ax = plt.subplots(figsize=(5,4))
  #plt.title('Number of refinements: '+str(count))
  ax.hist2d(D2BH, gas_r, bins=128, range=[[0, xaxis_max], [0,yaxis_max]])
  ax.plot(x, y, color='white')
  ax.axhline(min_size, 0, 1, color='yellow', linestyle='--')
  ax.axvline(1., 0, 1, color='yellow', linestyle='--')
  ax.set_xlim(0, xaxis_max)
  ax.set_ylim(0, yaxis_max)
  ax.set_xlabel('Distance from Center [kpc]')
  ax.set_ylabel('Cell Size Radius [kpc]')
  fig.savefig(output_name)
  plt.close()
