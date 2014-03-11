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

im = 100 # number of individuals immigrating from regional pool per time step

num_patches = 20 # number of patches in each local community

lgp = 0.99 # log-series parameter, typically approaches 1 for ecological
           # communities. represents underlying structure of regional pool.

state = 'heterogeneous' # homogeneity or heterogeneity among local communities
            # local communities vary along a single environmental axis, e.g.,
            # precipitation, exposure, runoff, etc.

time = 100 # number of time steps. 


northCOM, southCOM = cf.microbide(im, num_patches, lgp, state, time)

print len(northCOM),'patches in north and', len(southCOM),'patches in south'