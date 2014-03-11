from __future__ import division
import sys                                            
import numpy as np
import random
from random import randrange
import decimal
import cloud

sys.path.append("/Users/lisalocey/Desktop/evolution-canyon/microbide/models/coreFunctions")
import coreFunctions as cf


""" This script will run the microbide model and generate site by species
    matrices.
"""


get_siteBySpeciesMatrix(northCOM, southCOM):
    
    




num_patches = 20 # number of patches in each local community

lgp = 0.99 # log-series parameter, typically approaches 1 for ecological
           # communities. represents underlying structure of regional pool

state_list = ['heterogeneous','homogeneous'] # homogeneity or heterogeneity 
            # among local communities local communities vary along a single 
            # environmental axis, e.g.,mean daily temp, precipitation, etc.

time_list = [100,500,1000] # number of time steps 
im_list = [100, 200, 400]  # number of individuals immigrating
                           # from regional pool per time step

for state in state_list:
    for im in im_list:
        for time in time_list:
            
            northCOM, southCOM = cf.microbide(im, num_patches, lgp, state, time)
            print len(northCOM),'patches in north and',
            print len(southCOM),'patches in south'
            
            siteBySpeciesMatrix = get_siteBySpeciesMatrix(northCOM, southCOM)
            