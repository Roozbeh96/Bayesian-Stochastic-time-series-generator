# %% import packages
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import h5py
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

#%% Load Data
WT7VorX = sio.loadmat('G:/My Drive/Research/DATABASE/Vortex statistics/data/mesh_07ms_vortex_properties.mat')
WT7VorX_scales = sio.loadmat('G:/My Drive/Research/DATABASE/Vortex statistics/data/mesh_07ms_scales.mat')
VF_PIVWT7 = h5py.File(r'VF_PIVWT7(struc).mat', 'r')

#%%
VF_PIVWT7 = VF_PIVWT7['VF_PIVWT7'] 
VF_PIVWT7_x = np.array(VF_PIVWT7['x'])
VF_PIVWT7_z = np.array(VF_PIVWT7['z'])


# %%
for key in WT7VorX.keys():
    print(key)

zprograde = WT7VorX['z_pro']
zretrograde = WT7VorX['z_ret']



# %%
for key in WT7VorX_scales.keys():
    print(key)

zprof = WT7VorX_scales['z']
delta = WT7VorX_scales['delta']


# %%

numb_prograde = np.zeros((21,1))
numb_retrograde  = np.zeros((21,1))

for i in range(0,21):
    numb_prograde[i] = np.sum((zprograde > zprof[i]) & (zprograde < zprof[i+1]))
    numb_retrograde[i] = np.sum((zretrograde > zprof[i]) & (zretrograde < zprof[i+1]))


# %%
fig, axis = plt.subplots(1, 2, figsize=(10, 5))  

# Plotting for the first subplot
axis[0].plot(numb_prograde/(10000*np.max(VF_PIVWT7_x)**2*np.max(VF_PIVWT7_z)**2),

              zprof[0:21, 0] / delta[0, 0])
axis[0].set_xlabel('$\#$Number of prograde vortices/(Area)', fontsize=12)  
axis[0].set_ylabel(r'z/$\delta$', fontsize=12)  
axis[0].grid(True)  

# Plotting for the second subplot
axis[1].plot(numb_retrograde/(10000*np.max(VF_PIVWT7_x)**2*np.max(VF_PIVWT7_z)**2),
              zprof[0:21, 0] / delta[0, 0], 'r')
axis[1].set_xlabel(r'$\#$Number of retrograde vortices/(Area)', fontsize=12)  
axis[1].set_ylabel(r'z/$\delta$', fontsize=12)  
axis[1].grid(True)  

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()
# %%
