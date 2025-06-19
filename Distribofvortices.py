

'''Packages'''
#%%
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.io as sio
from sklearn.cluster import DBSCAN
from scipy.io import savemat
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

#%%
'''Data Reading'''

print("Start reading the Data...\n")
VF_PIVWT7 = h5py.File(r'VF_PIVWT7(struc).mat', 'r')
print("Data reading is completed.\n")

#%%
for key in VF_PIVWT7.keys():
    print(key)

VF_PIVWT7 = VF_PIVWT7['VF_PIVWT7'] 
for key in VF_PIVWT7.keys():
    print(key)

#%%
VF_PIVWT7_x = np.array(VF_PIVWT7['x'])
VF_PIVWT7_z = np.array(VF_PIVWT7['z'])
VF_PIVWT7_Lambda_ci = np.transpose(VF_PIVWT7['Lambda_ci'],(2,1,0))

#%%
# WT7VorX_scales = sio.loadmat('/Users/roozbehehsani/Google Drive/My Drive/Research/DATABASE/' 
                    # 'Vortex statistics/data/mesh_07ms_scales.mat')
WT7VorX_scales = sio.loadmat('G:/My Drive/Research/DATABASE/Vortex statistics/data/mesh_07ms_scales.mat')
for key in WT7VorX_scales.keys():
    print(key)

#%%    
lambda_prof = WT7VorX_scales['lambda']
z_prof = WT7VorX_scales['z']
delta = WT7VorX_scales['delta']
print(z_prof[0:22]/delta)
lambda_T_bar = np.mean(lambda_prof[0:22,0])
print(lambda_T_bar)

#%%
plt.figure
plt.plot(lambda_prof[0:22,0], z_prof[0:22]/delta)

# Add labels and title
plt.xlabel(r'$\lambda_{T}$')
plt.ylabel(r'z/$\delta$')
plt.show()

#%%
xint = int(np.floor(VF_PIVWT7_x[0,-1]/lambda_T_bar))

zint = int(np.floor(VF_PIVWT7_z[0,-1]/lambda_T_bar))

#%%
ix = 0
retro_vor_tot = []
prog_vor_tot = []
for S in range(VF_PIVWT7_Lambda_ci.shape[2]):
    for row in range(zint):
        # for col in range(xint):
        Matrix = VF_PIVWT7_Lambda_ci[row*(zint):(row+1)*(zint),
                            :,S]
        bin_Matrix = np.where(Matrix < 0, -1, np.where(Matrix > 0, 1, 0))
        clustered_matrix = np.zeros_like(bin_Matrix, dtype=int)
        
        pos_ind = np.argwhere(bin_Matrix == 1)
        if pos_ind.size > 0:
            db_1 = DBSCAN(eps=1, min_samples=1).fit(pos_ind)
            labels_1 = db_1.labels_
            retro_vor_tot.append([len(set(labels_1)) - (1 if -1 in labels_1 else 0),
                                np.mean(VF_PIVWT7_z[0,row*(zint):(row+1)*(zint)],axis = 0)])
            for idx, label in zip(pos_ind, labels_1):
                clustered_matrix[idx[0], idx[1]] = label + 1  # Start clusters at 1
        
        neg_ind = np.argwhere(bin_Matrix == -1)
        if neg_ind.size > 0:
            db_neg1 = DBSCAN(eps=1, min_samples=1).fit(neg_ind)
            labels_neg1 = db_neg1.labels_
            prog_vor_tot.append([len(set(labels_neg1)) - (1 if -1 in labels_neg1 else 0), 
                                np.mean(VF_PIVWT7_z[0,row*(zint):(row+1)*(zint)],axis = 0)])
            for idx, label in zip(neg_ind, labels_neg1):
                clustered_matrix[idx[0], idx[1]] = -(label + 1) 
        ix =+1
        # fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        # # Original matrix
        # axes[0].imshow(bin_Matrix, cmap='coolwarm', interpolation='none')
        # axes[0].set_title('Original Matrix')

        # # Clustered matrix
        # im2 = axes[1].imshow(clustered_matrix, cmap='tab20c', interpolation='none')
        # axes[1].set_title('Clustered Matrix')
        # fig.colorbar(im2, ax=axes[1])  # Add colorbar to the second plot
        # plt.show()
pos_vor_tot = np.array(retro_vor_tot)
neg_vor_tot = np.array(prog_vor_tot )

# %%
savemat('stat_vor_lambda_to_x_PIVWT7.mat', {'retrograde_vor_tot': retro_vor_tot, 'prograde_vor_tot': prog_vor_tot})


# %%
