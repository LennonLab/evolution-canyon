from __future__ import division 

import matplotlib # statements for running on IU servers be sure to comment out
matplotlib.use('Agg') # if running on local machine, resets where python looks

import matplotlib.pyplot as plt 
import sys, csv

import EClandscape as land
#import ECfigure as ECfig
import ECMicrobideCore as model
import ECfunctions as funx

#path = '~Desktop/evolution-canyon/microbide/SbyS/'
 

###########################  GET CONDITIONS  ################################### 

""" Code to runs the microbide model and generates site-by-species matrices."""

num_patches = 20 # number of patches on each side of Evolution Canyon (EC)
lgp = 0.92 # log-series parameter; underlying structure of regional pool

conditions = [['same', 'rand', 'rand'],  
             ['differ', 'rand', 'rand'],
             ['same',  'env',  'env'],
             ['differ', 'rand', 'env'],
             ['differ', 'env',  'rand'],
             ['differ',  'env', 'env']]
             
""" conditions is a list of modeling parameters for different conceptual 
    predictions representing extreme ends of a continuum of possible differences
    in environment and whether entering and exiting from dormancy has a
    stochastic component
    
    Conditions for...
        Conceptual prediction 1: 
            Environments have the same effect
            Exiting and entering dormancy has a large stochastic component 
    
        Conceptual prediction 2.
            Environments have different effects. 
            Exiting and entering dormancy has a large stochastic component
    
        Conceptual prediction 3.
            Environments have the same effect
            Exiting and entering dormancy is environmental
    
        Conceptual prediction 4.
            Environments have different effects. 
            Entering dormancy has a large stochastic component 
            Exiting dormnancy is environmental
    
        Conceptual prediction 5.
            Environments have different effects. 
            Entering is environmental
            Exiting dormancy has a large stochastic component
            
        Conceptual prediction 6.
            Environments have different effects. 
            Entering is environmental, Exiting is environmental
    
"""    
   
####################  GENERATE SITE BY SPECIES DATA  ###########################
for ic, combo in enumerate(conditions):
    
    envDiff, enterD, exitD = combo 

    landscapeLists = land.get_landscape(combo) # characterizing the landscape
    
    NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys = landscapeLists[0]
    SRow2Ys, Ncounts, Nverts, Scounts, Sverts = landscapeLists[1]
    
    COM = model.microbide(combo, Ncounts, Nverts, Scounts, Sverts, ic)
                            # run the model & return the communities
    
    nCOM, sCOM = funx.SpeciesInPatches(COM, NRowXs, NRow1Ys, NRow2Ys,
                                            SRowXs, SRow1Ys, SRow2Ys)                       
    
    
    SbyS = funx.get_SitebySpecies([nCOM, sCOM]) # get site by species matrix 
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
                                        
    path = '/N/dc2/projects/Lennon_Sequences/2014_EvolutionCanyon/microbide/SbyS/'
    #path = '~Desktop/evolution-canyon/microbide/SbyS/'
    
    fileName = 'Condition'+str(ic+1)
                                                                                
    fileName = 'Condition'+str(ic+1)
    OUT = open(path + fileName + '.txt','w')
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