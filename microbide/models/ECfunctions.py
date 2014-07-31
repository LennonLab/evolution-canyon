from __future__ import division 
import numpy as np
from scipy.stats import gaussian_kde
import sys


def get_kdens(summands): # Finds the kernel density function across a sample
    
    density = gaussian_kde(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)), len(summands))
    #xs = np.linspace(0,1, len(summands))
    density.covariance_factor = lambda : 0.1
    density._compute_covariance()
    
    return [xs, density(xs)]
    

def e_var(SAD):
    
    """ Calculates Smith and Wilson's evenness index.
    
    Smith, B. and J.B. Wilson (1996) A consumer's guide to evenness indices.
    OIKOS. 76, 70-82.
    """
    
    P = np.log(SAD)
    S = len(SAD)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/np.pi*np.arctan(X) 
    return(evar)
    
    
def get_match(Nverts, Sverts, x1, y1, Ncounts, Nmax, Scounts, Smax, opt1, opt2): 

    dmin = 10**6
    match = 0
    
    if x1 >= 0.5:
        
        for j, vert in enumerate(Nverts):
            x2, y2 = vert
                
            dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if dist < dmin:
                dmin = dist
                match = np.abs(opt1 - Ncounts[j]/Nmax) # mismatch to North niche
        
    elif x1 < 0.5:
        dmin = 10**6
        match = 0
        
        for j, vert in enumerate(Sverts):
            x2, y2 = vert
            
            dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if dist < dmin:
                dmin = dist
                match = np.abs(opt2 - Scounts[j]/Smax) # mismatch to South niche
    
    return match
    


def getNSDA(COM, verbose=True):
    
    """ A function to find total abundance, species richness, active abundance,
        dormant abundance, dominance, and evenness of a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    """
    
    Ss = []
    A = 0
    D = 0
    N = 0
    
    for patch in COM:
        N += len(patch)
        for i in patch:
            if i[1] > 0:
                D += 1
                Ss.append(i[0])
            elif i[1] == 0:
                A += 1
                Ss.append(i[0])
            
    uniqueSs = list(set(Ss))
    uniqueSs.sort()
    
    abVector = [] # vector of abundances ordered according to species identity
    AD = []
    for S in uniqueSs:
        ab = Ss.count(S)
        abVector.append([S, ab])
        AD.append(ab)
    S = len(list(set(Ss)))
    
    BP = max(AD)/N # Berger Parker index of dominance
    evar = e_var(AD) # Smith and Wilson's evenness index 
    
    if verbose == True:
        return [N, S, D, A, round(BP,3), round(evar,3), abVector]
    elif verbose == False:
        return abVector
        


def SpeciesInPatches(COM, NRowYs, NRow1Xs, NRow2Xs, SRowYs, SRow1Xs, SRow2Xs):
    
    NXs = NRow1Xs + NRow2Xs
    NYs = NRowYs + NRowYs
    
    SXs = SRow1Xs + SRow2Xs
    SYs = SRowYs + SRowYs
    
    NPatches = [list([]) for _ in range(len(NXs))]
    SPatches = [list([]) for _ in range(len(SXs))]
    
    for i, ind in enumerate(COM):
        x, y = ind[2], ind[3]
        
        for p, xcoord in enumerate(NXs):
            if xcoord - .02 <= x and xcoord + .02 >= x:         
                ycoord = NYs[p]
                if ycoord - .1 <= y and ycoord + .1 >= y:
                    NPatches[p].append([ind[0], ind[1]])
                    
            elif SXs[p] - .02 <= x and SXs[p] + .02 >= x:         
                ycoord = SYs[p]
                if ycoord - .1 <= y and ycoord + .1 >= y:
                    SPatches[p].append([ind[0], ind[1]])
                    
    return [NPatches, SPatches]           

  

def get_SitebySpecies(COMs):
    
    """ A function to generate the site by species matrix in a particular 
    formate required by Mario's PCoA analyses. 
    """
    
    nABV = getNSDA(COMs[0], False) # False directs the function to only return
    sABV = getNSDA(COMs[1], False) # the species abundance vector
    
    nmax = nABV[-1][0]
    smax = sABV[-1][0]
    maxID = max([nmax, smax]) # largest species id#
    SbyS = []
    
    for COM in COMs:
        for patch in COM:
            all_list = [0]*maxID
            active_list = [0]*maxID
            
            for ind in patch:
                
                spID = ind[0]
                all_list[spID-1]  += 1
                
                if ind[1] == 0:
                    active_list[spID-1] += 1
            
            SbyS.append(all_list) # DNA analogue
            SbyS.append(active_list) # RNA analogue
    
    # Remove any species with 0 abundance across sites
    SbyS = zip(*SbyS) # zip 
    SbyS = [list(row) for row in SbyS if any(x != 0 for x in row)] # remove 0's
    SbyS = zip(*SbyS) # un-zip
    
    S = len(SbyS)
    #print 'number sites:', len(SbyS)
    
    if len(SbyS)%2 > 0:
        print 'unequal number of sites between Active and All'
        sys.exit()
    
    siteNumbers = []
    for i in range(int(S/2)):
        siteNumbers.extend([i+1, i+1])
    
    minRowN = 10**8
    minRow = list()
    
    for i, row in enumerate(SbyS):
        row = list(row)
        name = str()
            
        if i <= 39: name = 'north'
        else: name = 'south'
        
        if (i+1)%2 > 0:
            SbyS[i] = [0.03, name + str(siteNumbers[i]) + '_all', S] + row
        else:            
            SbyS[i] = [0.03, name + str(siteNumbers[i]) + '_active', S] + row
            
        if sum(row) <= minRowN:
            minRowN = sum(row)
            minRow = SbyS[i]
        
    
    print 'smallest site abundance:', minRowN,', site name:', minRow[1] 
    return SbyS  
