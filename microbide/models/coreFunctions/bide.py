from __future__ import division
import sys                                            
import numpy as np
from random import randrange, choice
import DiversityTools as dt


""" Core functions for Evolution Canyon Project modeling """

def get_param(ID, param_dict):
    """ A function to set species specific parameter values """
    param = np.random.uniform(0.01, 0.99)
    param_dict[ID] = param
        
    return param_dict


def reproduction(nCOM, sCOM):
    
    """ A function to simulate reproduction in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    for COM in [nCOM, sCOM]:
        for i, patch in enumerate(COM):
            while patch:    
                ind = randrange(len(patch))
                # choose an individual at random from patch p 
                if patch[ind][1] == 0: # keep going until an active one is found
                    COM[i].append(patch[ind])
                    break
                                                                
    return [nCOM, sCOM]



def dormancy(nCOM, sCOM, enVal, optima_dict):
    
    """ A function to simulate the transition to or from dormancy. Here, an
    imperfect match between the environment and the species optimum can still
    result in the transition from dormancy to activity. A near perfect match
    between the environment and the species optimum can also result in the 
    transition from activity to dormancy (but is very unlikely).
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    for COM in [nCOM, sCOM]:
        for i, patch in enumerate(COM):
            p = randrange(len(COM)) # choose a patch from COM at random and   
            j = randrange(len(COM[p])) # can individual at random from patch p
            ind = COM[p][j]
        
            optima = optima_dict[ind[0]] # the species environmental optima 
            Err = np.abs(optima - enVal) # difference between species
                                    # optima and patch environment;   
                                    # value of 0 means no disagreement.
                                        
            if np.random.binomial(1, Err) == 0: COM[p][i][1] = 0 
            # higher error means higher chance of being dormant
            else: COM[p][i][1] += 1
     
    return [nCOM, sCOM]
      
          

def death_emigration(nCOM, sCOM, enVal, optima_dict, ceiling=50):

    """ A function to simulate death/emigration in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    for COM in [nCOM, sCOM]:
        for index, patch in enumerate(COM):
            
            while len(patch) > ceiling:    
            
                i = randrange(len(patch)) # choose individual at random from patch p
                ind = patch[i]
    
                optima = optima_dict[ind[0]] # the species environmental optima 
                Err = np.abs(optima - enVal) # difference between species
                                    # optima and patch environment;   
                                    # value of 0 means no disagreement.
                                    
                if np.random.binomial(1, Err) == 1: # higher error means
                    patch.pop(i)                # higher chance of leaving/death
                                        
        COM[index] = patch
    return [nCOM, sCOM]



def disperse(nCOM, sCOM, enVal, disp = 0.1):

    """ A function to simulate disperal between two local communities/areas.
    
    COM1 & COM2  :  the local communities; vectors of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    """
    
    for COM in [nCOM, sCOM]:
        for index, patch in enumerate(COM):
            p = randrange(len(COM))  # choose a patch at random from COM1
            i = randrange(len(COM[p])) # choose an individual at random
            ind = COM[p].pop(i)
            
            i = randrange(len(COM))               # Choose an individual to disperse                                 
            x = np.random.normal(loc=COM[i][2], scale=disp)     
	    y = np.random.normal(loc=COM[i][3], scale=disp)
	
	    if x > 1: x = 1 # keeping the individual within bounds
	    elif x < 0: x = 0
	    if y > 1: y = 1
	    elif y < 1: y = 1
	    COM[i][2], COM[i][3] = [x, y]
	        
            
            p = randrange(len(COM)) # choose a patch at random from COM2
            COM[p].append(ind)
        
    return nCOM, sCOM




def immigration(COMs, optima_dict, env_optima, decay_dict, dorm_decay, im):
    
    for i, COM in enumerate(COMs):
        for j, site in enumerate(COM):     
            propagules = np.random.logseries(0.9, im)  # list of propagules 
    
            for prop in propagules:
                if prop not in optima_dict:
                    optima_dict = get_param(prop, optima_dict)
                    decay_dict = get_param(prop, decay_dict)
                
                COMs[i][j].append([prop, choice([0, 1])]) # 0 is active, 1 is dormant
                                               
    return COMs