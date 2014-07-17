from __future__ import division
import sys, csv
import numpy as np
import random
from random import randrange
import decimal
import cloud

sys.path.append("/Users/lisalocey/Desktop/evolution-canyon/microbide/models/coreFunctions")
import coreFunctions as cf


""" Code to runs the microbide model and generates site-by-species matrices."""

num_patches = 20 # number of patches on each side of Evolution Canyon (EC)
lgp = 0.92 # log-series parameter; underlying structure of regional pool

""" conditions is a list of modeling parameters for different conceptual 
    predictions representing extreme ends of a continuum of possible differences
    in environment, dormancy, and secondary model parameters (e.g. dispersal,
    scale, etc.)
    
    condition 1. Does the environment differ between N and S? If yes, envs. will
    be polar opposites (0 and 1.0)
     
    condition 2. Does dormancy occur? 
    
    
    
    condition 3. Does dormancy differ 
    
    
    Conceptual prediction 1. 
    Conditions: Environments don't differ. Dormancy is random switching. 
    
    Conceptual prediction 2.
    Conditions: Environments differ. Dormancy is random switching.
    
    Conceptual prediction 3.
    Conditions: Environment do not differ. Dormancy is environmentally driven.
    
    Conceptual predictions 4.
    Conditions: Environments differ. Entering dormancy is random. Exiting dormancy
    is environmental. Classic bet-hedging 
    
    Conceptual predictions 5.
    Conditions. Environments differ. Entering is environmental. Existing is random...Scout hypothesis
    
    Conceptual predictions 6.
    Conditions: Environment differ. 
    
    """


conditions = []

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