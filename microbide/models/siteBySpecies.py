from __future__ import division
import sys, csv
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

num_patches = 20 # number of patches in each local community

lgp = 0.99 # log-series parameter, typically approaches 1 for ecological
           # communities. represents underlying structure of regional pool

state_list = ['heterogeneous']#,'homogeneous'] # homogeneity or heterogeneity 
            # among local communities local communities vary along a single 
            # environmental axis, e.g.,mean daily temp, precipitation, etc.

im_list = [300]  # number of individuals immigrating
                           # from regional pool per time step

time = 70
for state in state_list:
    for im in im_list:
        
        northCOM, southCOM = cf.microbide(im, num_patches, lgp, state, time)
        print len(northCOM),'patches in north and',
        print len(southCOM),'patches in south'
            
        SbyS = cf.get_SitebySpecies([northCOM, southCOM])
        S = len(SbyS[0]) - 3
        
        path = '/Users/lisalocey/Desktop/evolution-canyon/microbide/models'
        path = path + '/SiteBySpecies/'
        fileName = 'SiteBySpecies_' + str(state) + 'immigration=' + str(im)
        
        OUT = open(path + fileName + '.share','w+')
        writer = csv.writer(OUT, delimiter='\t')
        
        linedata = ['label', 'Group', 'numOtus']
        for i in range(S):
            linedata.append('Otu'+str(i))
        writer.writerow(linedata)
                
        for row in SbyS:
            writer.writerow(row)
    
        sys.exit()
        
            
            