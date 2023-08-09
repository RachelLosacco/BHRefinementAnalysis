# BH accretion rate (BHAR) vs. Star formation rate (SFR)
# Using DataLoader (python 3)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader

print('Starting')
# Specifications
version = 'debugging_level7'
path = '/orange/paul.torrey/rlosacco/'+version+'/'
enclosed_radius = [50., 1., 0.1, 0.01] #kpc
mycolors = ['red', 'green', 'blue', 'purple']
labels = ['<50kpc', '<1kpc', '<100pc', '<10pc']
snapnum = 171
a = 0.99881612
box = 20.  # Mpc^3
h = 0.6909 # from param.txt

# Create plot
fig, ax = plt.subplots(figsize=(5.5,4.5))
logylimits = [-3., 2.]
logxlimits = [-3., 2.]
xlimits = [1e-3, 1e2]
ylimits = xlimits
width = 1

# Constants
z = (1./a) - 1.
volume = (box/h)**3.
sim_mass_unit = 1e10/h  # simulation mass units


# Read in data
for i in range(len(enclosed_radius)):
  print('enclosed radius of ', labels[i])
  total_sfr = []
  total_bhar = []
  count = 0
  for groupn in range(50):
    print('group no. ', groupn)
    BHcat = DataLoader(path, part_types=[5], snap_num=snapnum, keys=['BH_Mdot', 'Coordinates'], sub_idx=0, fof_idx=groupn)
    BHAR = BHcat['BH_Mdot'] * 1e10 / 1e9   # original units: (1e10Msun/h) / (0.978Gyr/h)
    
    gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=['StarFormationRate', 'Coordinates'], sub_idx=0, fof_idx=groupn)
    SFR = gascat['StarFormationRate']
    
    # Calculate distance to nearest blackhole by hand
    gaspos = gascat['Coordinates']
    bhpos = BHcat['Coordinates']
    if len(bhpos) == 0:
      print('no bh')   # not all snapshots have a blackhole
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
      gasdist = np.sqrt(dx*dx + dy*dy + dz*dz)
    
    # Total SFR and BHAR within enclosed radius
    sumofsfr = np.sum(SFR[gasdist < enclosed_radius[i]])
    total_sfr.append(sumofsfr)
    total_bhar.append(np.sum(BHAR))
      
    if sumofsfr > 0:
      count += 1
    
  #Plotting
  print('len of sfr', len(total_sfr))
  print('len of bhar', len(total_bhar))
  print('count=', count)
  ax.scatter(total_sfr, total_bhar, s=2, color=mycolors[i], label=labels[i])
ax.set_xlabel("Star Formation Rate [$M_\odot$/yr]")
ax.set_ylabel("Black Hole Accretion Rate [$M_\odot$/yr]")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(ylimits)
ax.set_xlim(xlimits)
x = np.linspace(xlimits, 100)
ax.plot(x, x, linestyle='solid', linewidth=width, color='k', alpha=.5)
ax.plot(x, 0.1*x, linestyle='dashdot', linewidth=width, color='k', alpha=.5)
ax.plot(x, 0.01*x, linestyle='dashed', linewidth=width, color='k', alpha=.5)
ax.plot(x, 0.001*x, linestyle='dotted', linewidth=width, color='k', alpha=.5)
ax.legend()
fig.savefig('BHARvSFR_snap'+str(snapnum)+version+'.png', bbox_inches='tight')
plt.close()

print('Done')

