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

sp_env_optima = 'random'
sp_dorm_decay = 'random'
northVal = 0.2
southVal = 0.8 
btn = True # dispersal between patches

northCOM, southCOM = cf.microbide(num_patches, lgp, northVal, southVal, sp_env_optima, sp_dorm_decay, btn)                                        
SbyS = cf.get_SitebySpecies([northCOM, southCOM]) 
S = len(SbyS[0]) - 3
 
                     
r1 = len(SbyS[0])
for i, row in enumerate(SbyS):
    r2 = len(row)
    if i%2 > 0 and sum(row[2:]) == 0:
        print 'there are no individuals in row',i
    if r1 != r2:
        print 'unequal sized rows in Site by Species matrix'
        sys.exit()
    r1 = r2
        
                                
path = '/Users/lisalocey/Desktop/evolution-canyon/microbide/models/SbyS/'
fileName = 'Cartoon1'

OUT = open(path + fileName + '.share','w')
writer = csv.writer(OUT, delimiter='\t')
                    
linedata = ['label', 'Group', 'numOtus']
for i in range(S):
    linedata.append('Otu'+str(i))
                    
writer.writerow(linedata)
                
for row in SbyS:
    if len(row) != r1:
        print 'row length has been corrupted'
        sys.exit()
                    
    writer.writerow(row)
    #print row
                        
OUT.close()