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

conditions = [[False, False, False],  
             [True, False, False],
             [False, True, True],
             [True, False, True],
             [True, True, False]]
             
""" conditions is a list of modeling parameters for different conceptual 
    predictions representing extreme ends of a continuum of possible differences
    in environment and whether entering and exiting from dormancy has a
    stochastic component
    
    Conditions for:
        Conceptual prediction 1: 
            Environments do not differ
            Dormancy is random switching
    
        Conceptual prediction 2.
            Environments differ. 
            Dormancy is random switching
    
        Conceptual prediction 3.
            Environments do not differ
            Dormancy is environmentally driven
    
        Conceptual prediction 4.
            Environments differ. 
            Entering dormancy has a stochastic component 
            Exiting is environmental (seed setting)
    
        Conceptual prediction 5.
            Environments differ
            Entering is environmental
            Exiting dormancy has a stochastic component
    
"""

for combo in conditions:
    envDiff, enterD, exitD = combo 
    
    nVal = 0.1 # environment of the north face
    sVal = 0.9 # environment of the south face

    nCOM, sCOM = cf.microbide(num_patches, lgp, nVal, sVal, 
                envDiff, enterD, exitD) # run the model & return the communities
                                                                                                                      
    SbyS = cf.get_SitebySpecies([nCOM, sCOM]) # get site by species matrix 
    S = len(SbyS[0]) - 3 # first 3 columns are non-species data
                 
                 
    r1 = len(SbyS[0])
    for i, row in enumerate(SbyS): # checking to ensure correct format for SbyS
        r2 = len(row)
        
        if i%2 > 0 and sum(row[2:]) == 0: # first 3 columns are non-species data
            print 'there are no individuals in row', i
        
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
                        
    OUT.close()