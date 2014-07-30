from __future__ import division 
from random import choice, randrange
import numpy as np

import ECfunctions as funx



def immigration(COM, oDict, dDict, im):    
    
    propagules = np.random.logseries(0.97, im)  # list of propagules 
    for p in propagules:
        if p not in oDict:
            oDict[p] = np.random.uniform(0.001, 0.99, 1).tolist()[0]
            dDict[p] = np.random.uniform(0.111, 0.99, 1).tolist()[0]
            
        x = np.random.uniform(0.001, 0.99, 1).tolist()[0]
        y = np.random.uniform(0.001, 0.99, 1).tolist()[0]
        
        COM.append([p, choice([0, 1]), x, y]) # 0 is active, 1 is dormant
                                                                                                                                 
    return COM


######################### COMMUNITY SIMULATION FUNCTION ########################
def microbide(combo, Ncounts, Nverts, Scounts, Sverts, ic):
    
    envDiff, enterD, exitD = combo 
    
    Nmax = max(Ncounts) # max heat value
    Smax = max(Scounts) # max heat value
    
    oDict, dDict, COM = [ {}, {}, [] ]
    
    N = 2*10**5 # Starting total abundance across the landscape
    COM = immigration(COM, oDict, dDict, N)
    
    for t in range(2*10**6):
        
        """ Immigration """
        x = np.random.binomial(1, 0.2)
        if x == 1:
            COM = immigration(COM, oDict, dDict, 1)
    
    
        """ Death """
        i = randrange(len(COM))
        spID, state, x1, y1, opt = [COM[i][0], COM[i][1], COM[i][2], 
                                            COM[i][3], oDict[COM[i][0]]]
        
        match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                                    Nmax, Scounts, Smax, opt)
        
        if state == 0:
            x = np.random.binomial(1, match)
            if x == 0: COM.pop(i) # likely to die from a low match    
	
	elif state == 1:
            d = dDict[spID]
            x = np.random.binomial(1, match) # multiplicative interaction between dormancy decay and mismatch to env.
            if x == 0: 
                x = np.random.binomial(1, d) # 'second chance' from dormancy
                if x == 1:
                    COM.pop(i) # die in dormancy, inherently neg. binomial  
           		
	""" Reproduction """
	i = randrange(len(COM))
        ind = COM[i]
        spID, state, x1, y1, opt = [ind[0], ind[1], ind[2], 
                                            ind[3], oDict[ind[0]]]
        
        match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                                    Nmax, Scounts, Smax, opt)
        
        x = np.random.binomial(1, match)
        if state==0 and x==1: COM.append(COM[i]) # reproduce from a high match   
	
	    
	""" Transition between activity and dormancy """
	i = randrange(len(COM))
        ind = COM[i]
        spID, state, x1, y1, opt = [ind[0], ind[1], ind[2], 
                                            ind[3], oDict[ind[0]]]
        
        if state == 0: # if active
            match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts, 
                                    Nmax, Scounts, Smax, opt)
            
            if enterD == 'random':
                match1 = 0
	        while match1 <=0 or match1 >=1: 
	           match1 = np.random.normal(match, 0.1, 1)
	        match = float(match1)
	    
	    x = np.random.binomial(1, match)
	    if x == 0: COM[i][1] = 1 # go dormant with a low match              
	
	if state == 1: # if dormant
	    match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts, 
                                    Nmax, Scounts, Smax, opt)
            
            if exitD == 'random': 
	        match1 = 0
	        while match1 <=0 or match1 >=1: 
	            match1 = np.random.normal(match, 0.1, 1)
	        match = float(match1)
	    
	    x = np.random.binomial(1, match)
	    if x == 1: COM[i][1] = 0 # activate with a high match   
	    
	                          
	""" Dispersal """
        i = randrange(len(COM))
        x = np.random.normal(loc=COM[i][2], scale=0.1)     
	y = np.random.normal(loc=COM[i][3], scale=0.1)
	
	if x > 1 or x < 0: COM.pop(i)
	elif y > 1 or y < 1: COM.pop(i)
	 
        t += 1
        if t%10000 == 0:
            print 'condition',ic+1,' time:',t, 'N:', len(COM)            
    
    return COM    
