# Fourier transform of gas inflow across a constant radius
# A0 and A2 are harmonics of the fourier
# See Beane et al. 2022, Fig 5
# See fourier_time.py and fourier_radius.py for more analysis
# Using DataLoader (python 3)

# Imports
import numpy as np
import matplotlib.pyplot as plt
from readData.DataLoader import DataLoader
import sys
from analysis import analyze
from copy import copy
import visualization.contour_makepic as cmakepic
import util.calc_hsml as calc_hsml
from scipy.interpolate import interp2d
import astropy.units as u
import astropy.constants as c


# Simulation Info
runtype = 'Spatial_Refine'
version = 'v0_cont'
snapnum = 68
size_min = 300

# Other info
constant_radius = 10     # kpc, fixed radius
constant_snap = 10       # z=0, fixed time
number_points = 200      # number of points around the circle
pixels = 1024
fov = 7.2                # kpc, fov-radius
bhrange = range(50, 101) # looping over these indeces

# Functions
## Mass Map
def get_massmap(cat, pixels, fov=30, face_on=True, edge_on=False, part_types=[4], mini=False):
    '''
    Input: Gas catalog from AREPO
    Returns: 2D mass distribution and image using cmakepic.simple_makepic
    '''
    if type(fov) == type([]):
        xrange=[-fov[0],fov[0]]
        yrange=[-fov[1],fov[1]]
        maxfov = np.max(fov)
    else:
        xrange=[-fov,fov]
        yrange=[-fov,fov]
        maxfov = fov

    coords = copy(cat['Coordinates'])/cat.h*cat.time
    vels = copy(cat['Velocities'])      # km/s
    masses = cat['Masses']*1e10/cat.h   # Msol
    gal_pos = cat['GroupPos']/cat.h*cat.time
    gal_vel = cat['GroupVel']
    
    # mini vs. full snapshots for smoothing length
    if mini:
      hsml = calc_hsml.get_particle_hsml(cat['Coordinates'][:,0], cat['Coordinates'][:,1], cat['Coordinates'][:,2])
    else:
      hsml = cat['SubfindHsml']

    box_cut = analyze.get_box_cut(coords, gal_pos, maxfov)
    coords -= gal_pos
    vels -= gal_vel
    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on)

    massmap, image = cmakepic.simple_makepic(coords[box_cut,0], coords[box_cut,1], weights=masses[box_cut], hsml=hsml[box_cut], xrange=xrange, yrange=yrange, pixels=pixels, set_maxden = 5.0e7, set_dynrng = 1.0e4)

    return massmap, image
    
## Momentum density map
def get_mommap(cat, pixels, bhpos, fov=30, face_on=True, edge_on=False, part_types=[4], mini=False):
    '''
    Input: Gas catalog from AREPO
    Returns: 2D momentum distribution and image using cmakepic.simple_makepic
    '''
    if type(fov) == type([]):
        xrange=[-fov[0],fov[0]]
        yrange=[-fov[1],fov[1]]
        maxfov = np.max(fov)
    else:
        xrange=[-fov,fov]
        yrange=[-fov,fov]
        maxfov = fov

    coords = copy(cat['Coordinates'])/cat.h*cat.time
    vels = copy(cat['Velocities'])              # km/s
    masses = cat['Masses']*1e10/cat.h           # Msol
    gal_pos = cat['GroupPos']/cat.h*cat.time
    gal_vel = cat['GroupVel']
    
    # r-hat and velocities to make radial velocity
    gasx = coords[:,0]
    gasy = coords[:,1]
    gasz = coords[:,2]
    bhx = bhpos[0,0]
    bhy = bhpos[0,1]
    bhz = bhpos[0,2]
    dx = gasx - bhx
    dy = gasy - bhy
    dz = gasz - bhz
    dist = np.sqrt(dx*dx + dy*dy + dz*dz)
    radial_unit_vector = [dx/dist, dy/dist, dz/dist]
    radial_unit_vector = np.asarray(radial_unit_vector)
    radial_unit_vector = np.transpose(radial_unit_vector)
    rxvx = radial_unit_vector[:,0] * vels[:,0]
    ryvy = radial_unit_vector[:,1] * vels[:,1]
    rzvz = radial_unit_vector[:,2] * vels[:,2]
    vel_r = rxvx + ryvy + rzvz    # km/s
    
    # mini vs. full snapshots for smoothing length
    if mini:
      hsml = calc_hsml.get_particle_hsml(cat['Coordinates'][:,0], cat['Coordinates'][:,1], cat['Coordinates'][:,2])
    else:
      hsml = cat['SubfindHsml']

    box_cut = analyze.get_box_cut(coords, gal_pos, maxfov)
    coords -= gal_pos
    vels -= gal_vel
    coords, vels = analyze.get_rotate_data(coords, vels, masses, face_on=face_on, edge_on=edge_on)
    
    # weighting by momentum gives momentum density map [Msol/km/s]
    mommap, image = cmakepic.simple_makepic(coords[box_cut,0], coords[box_cut,1], weights=masses[box_cut]*vel_r[box_cut], hsml=hsml[box_cut], xrange=xrange, yrange=yrange, pixels=pixels, set_maxden = 5.0e7, set_dynrng = 1.0e4)
    return mommap, image

## Calculate information
def get_info(runtype, version, snapnum, n, pixels, fov, number_points, constant_radius):
    '''
    Input: Information about the path, bh index in bhrange, snapshot number, and default values
    Output: if there is a bh or not (boolean), A0, A2, A2/A0 ratio, gas inflow value, and scale factor of snapshot
    '''
    # Get catalogues 
    path = '/orange/paul.torrey/rlosacco/'+runtype+'/'+version+'/'
    outname = '/home/rlosacco/'+runtype+'/Plotting/fourier_gas_'+version+'_snap'+str(snapnum).zfill(3)+'_fov'+str(fov)+'_'+str(n)+'.png'
    gaskeys = ["Coordinates", "Velocities", "Masses", "GroupPos", "GroupVel", "SubfindHsml"]
    bhkeys = ["Coordinates"]
    
    gascat = DataLoader(path, part_types=[0], snap_num=snapnum, keys=gaskeys, sub_idx=0, fof_idx=n)
    bhcat = DataLoader(path, part_types=[5], snap_num=snapnum, keys=bhkeys, sub_idx=0, fof_idx=n)
    
    gaspos = gascat['Coordinates']/gascat.h*gascat.time
    bhpos = bhcat['Coordinates']/bhcat.h*bhcat.time
    gasvel = gascat['Velocities']
    gasmass = gascat['Masses']*1e10/gascat.h  # Msol
    galpos = gascat['GroupPos']/gascat.h*gascat.time
    
    # is there a bh?
    if len(bhpos) == 0:
        print('no bh')
        nobh = True
        a0 = 1
        a2 = 0 
        mass_flow = 0
        #sys.exit()
    else:    
        # Distances and angles
        nobh = False
        gasx = gaspos[:,0]
        gasy = gaspos[:,1]
        gasz = gaspos[:,2]
        bhx = bhpos[0,0]
        bhy = bhpos[0,1]
        bhz = bhpos[0,2]
        dx = gasx - bhx
        dy = gasy - bhy
        dz = gasz - bhz
        dist = np.sqrt(dx*dx + dy*dy + dz*dz)  # kpc?
        phi = np.arctan2(dy, dx)
    
        # Gas map
        if "SubfindHsml" not in gascat:
            print("mini")
            gas_map, gas_image = get_massmap(gascat, pixels, fov=fov, part_types=0, mini=True)        # Msol
            mom_map, mom_image = get_mommap(gascat, pixels, bhpos, fov=fov, part_types=0, mini=True)  # Msol/km/s
        else:
            gas_map, gas_image = get_massmap(gascat, pixels, fov=fov, part_types=0)
            mom_map, mom_image = get_mommap(gascat, pixels, bhpos, fov=fov, part_types=0)
    
        vel_map = mom_map/gas_map   # 2D velocity map [km/s]
        
        # Function for surface density
        ## maybe fix linspace to hit the center of the pixels instead of the edges
        dist_per_pixel = np.linspace(-fov, fov, pixels, endpoint=True)       # kpc
        surface_density = interp2d(dist_per_pixel, dist_per_pixel, gas_map)  # Msol/kpc2
        velocity_distr = interp2d(dist_per_pixel, dist_per_pixel, vel_map)   # radial velocity density, km/s/kpc2
        
        # Calculating A0 and A2
        a0_real = 0.
        a2_real = 0.
        a2_img = 0.
        mass_flow = 0.
        for i in range(number_points):
            phi_i = i*2.*np.pi/number_points
            x_i = constant_radius * np.cos(phi_i)  # convert from polar to cartesian coordinates
            y_i = constant_radius * np.sin(phi_i)  # kpc
            rho_i = surface_density(x_i, y_i)[0]   # Msol/kpc2
            # A0 and A2  
            a0_real += rho_i
            a2_real += rho_i * np.cos(2.*phi_i)
            a2_img += rho_i * np.sin(2.*phi_i)
            # dM/dt (inflow)
            vel_i = velocity_distr(x_i, y_i)[0]                      # km/s/kpc2
            dMdt = rho_i * (u.M_sun / u.kpc**2.) * vel_i * (u.km / u.s / u.kpc**2.) * constant_radius * (u.kpc) * 2*np.pi   # (Msol/kpc2)(km/s/kpc2)(kpc)
            mass_flow += dMdt.decompose()
        
        a0 = a0_real
        a2 = np.sqrt(a2_real**2. + a2_img**2.)
    
    ## Plot
    print('plotting, saving as', outname)
    deltax = (galpos[0] - bhx)
    deltay = (galpos[1] - bhy)
    print('delta bhpos', deltax, deltay)
    circle = plt.Circle((deltax,deltay), constant_radius, color='r', fill=False)
    fig, ax = plt.subplots(figsize = (15.0, 15.0))
    ax.imshow(gas_image, origin='lower', extent=(-fov, fov, -fov, fov)) #, vmax=10, vmin=5)
    ax.plot
    ax.add_patch(circle)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.subplots_adjust(bottom=0, top=1, left=0, right=1)
    fig.savefig(outname, dpi=pixels/2. )
    plt.close()
    
    return nobh, a0, a2, a2/a0, mass_flow, gascat.time

a0s = []
a2s = []
a2a0s = []
mdots = []
scalefactors = []
for n in bhrange:
    print('galaxy', n)
    nobh, a0, a2, a2a0, mass_inflow, scalefactor = get_info(runtype, version, snapnum, n, pixels, fov, number_points, constant_radius)
    if nobh == True:
        print('no bh')
        continue
    else:
        a0s.append(a0)
        a2s.append(a2)
        a2a0s.append(a2a0)
        mdots.append(mass_inflow.value)
        scalefactors.append(scalefactor)
            
# Plotting
fig, ax = plt.subplots(figsize=(5, 4))
ax.hist(a2a0s, bins=50, range=(min(a2a0s), max(a2a0s)), histtype='step')
ax.set_xlabel("|A2/A0|")
fig.savefig('fourier_a2a0shist_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(5, 4))
ax.hist(mdots, bins=50, range=(min(mdots), max(mdots)), histtype='step')  #.to(u.M_sun/u.kpc/u.km/u.s)
ax.set_xlabel("dM/dt [kg/(m2 s)]")
fig.savefig('fourier_mdotshist_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(a2a0s, mdots)
ax.set_xlabel("|A2/A0|")
ax.set_ylabel("dM/dt [kg/(m2 s)]")
fig.savefig('fourier_a2a0smdots_'+version+'_snap'+str(snapnum)+'.png', bbox_inches='tight')
plt.close()


with open('fourier_info_'+version+'_snap'+str(snapnum)+'.txt', 'w') as outfile:
    for a0, a2, a2a0, mdot in zip(a0s, a2s, a2a0s, mdots):
        outfile.write(f'{a0}\t {a2}\t {a2a0}\t {mdot}\n')

