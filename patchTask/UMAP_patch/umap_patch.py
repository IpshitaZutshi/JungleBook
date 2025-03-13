# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 12:37:51 2024

@author: Laura Ribalta
"""
# from utils import *
import scipy
from scipy import io
import pandas as pd 
import csv
import math
import glob
import mat73
import os
import random
import umap
from umap import UMAP
#import matplotlib as plt

basepath = r'Z:\Buzsakilabspace\LabShare\ZutshiI\patchTask\N7\N7_241205_sess17'

spikes_file = glob.glob(basepath + '/*.data.mat')
behav_file = glob.glob(basepath + '/*.position_behavior_speed.mat')

matData = scipy.io.loadmat(spikes_file[0])

spikes_data = matData['data']

spikes_data.shape



random.seed(10)
#umap.umap_.UMAP
##################  Perform UMAP dimensionality reduction ########################
utrans = UMAP(n_neighbors=20, n_components=6, metric='cosine', metric_kwds=None, output_metric='euclidean', 
                output_metric_kwds=None, n_epochs=None, learning_rate=1.0, init='spectral', min_dist=0.1, spread=1.0, low_memory=True, 
                n_jobs=-1, set_op_mix_ratio=1.0, local_connectivity=1.0, repulsion_strength=1.0, negative_sample_rate=5, transform_queue_size=4.0, 
                a=None, b=None, random_state=None, angular_rp_forest=False, target_n_neighbors=-1, target_metric='categorical', target_metric_kwds=None, target_weight=0.5, 
                transform_seed=42, transform_mode='embedding', force_approximation_algorithm=False, verbose=False, tqdm_kwds=None, unique=False, densmap=False, dens_lambda=2.0, 
                dens_frac=0.3, dens_var_shift=0.1, output_dens=False, disconnection_distance=None, precomputed_knn=(None, None, None))
#movetimes1 = np.arange(0,len(sspikes1[:,0]),5)
#print(len(movetimes1))

uproj1 = utrans.fit_transform(spikes_data[:,:])
## fit transform behavior only

#import matplotlib.pyplot as plt

# plt.hsv()
# for deg1 in [-80, -20, 20,80,50,-50]:
#     for deg2 in [20]:
#         fig = plt.figure(figsize= (20,20), dpi = 300)

#         ax = fig.add_subplot(111, projection = '3d')       
#         ax.axis('off')
#         ax.set_title(str(deg1) + '_' + str(deg2))
#         ax.scatter(uproj1[:,0], 
#                    uproj1[:,1],
#                    uproj1[:,2],c='Crimson')

#                    #c = xx_tmaze[speed_tmaze>sp][movetimes1],
#                   #s = 10)
#         ax.view_init(deg1,deg2)
#         plt.show()



df = pd.DataFrame(uproj1)

# saving the dataframe
name_file = behav_file[0].split('\\')[-1]
name_file = os.path.splitext(name_file)[0]

df.to_csv(basepath + os.sep+ 'Umap_'+ name_file+'.csv',header=False, index=False)
#         pickle.dump(Umap_kNN, open(basepath + os.sep + basename + "."+ load_name  +"_"+ str(n_component)  +".p", "wb" ))

