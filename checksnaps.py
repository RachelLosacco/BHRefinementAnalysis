# Check to see if the outputted snapshot scale factors 
# line up with the desired scale factors
# My snapshots were being skipped over because the time steps
# were smaller than the snapshot frequency; this was to 
# check while debugging that issue

# Imports
import numpy as np
import matplotlib.pyplot as plt
import h5py

# Directory path and snapshots
version = 'vanilla6'
snapnum = 160
output_list_filename = '/home/rlosacco/vanilla/eighthundredsnaps.txt'
orange_path = '/orange/paul.torrey/rlosacco/'+version+'/'

# Get scale factors from output list (desired)
desired_a = np.genfromtxt(output_list_filename)[:,0]

# Get scale factors from snapshots (actual)
actual_a = []
for snap in range(snapnum):
    file_name = orange_path+'snapdir_'+str(snap).zfill(3)+'/snapshot_'+str(snap).zfill(3)+'.0.hdf5'
    with h5py.File(file_name, "r") as ofile:
        header = ofile['Header']
        a = header.attrs['Time']
        actual_a.append(a)
        
print('actual', len(actual_a))
print('desired', len(desired_a))

print('last few times of actual: ', actual_a[-3:])

# Plotting
fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(range(len(actual_a))[4:], actual_a[4:], label='Actual')
ax.plot(range(len(desired_a))[4:], desired_a[4:], label='Desired')
ax.legend()
ax.set_xlabel('Snapshot Number')
ax.set_ylabel('Scalefactor at Snapshot')
fig.savefig('vanilla_checksnaps.png')
plt.close()
