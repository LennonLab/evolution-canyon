from __future__ import division
import sys                                            
import numpy as np
import random
from random import randrange
import decimal



def get_SitebySpecies(COMs):
    
    """ A function to generate the site by species matrix in a particular 
    formate required by Mario's PCoA analyses. 
    """
    
    nABV = getNSDA(COMs[0], False) # False directs the function to only return
    sABV = getNSDA(COMs[1], False) # the species Id-abundance vector
    
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
                
                if ind[1] == 'a':
                    active_list[spID-1] += 1
            
            SbyS.append(all_list) # DNA analogy
            SbyS.append(active_list) # RNA analogy
    
    SbyS = zip(*SbyS) # these 
    SbyS = [list(row) for row in SbyS if any(x != 0 for x in row)]
    SbyS = zip(*SbyS)
    
    S = len(SbyS[0])
    if len(SbyS)%2 > 0:
        print 'unequal number of patches'
        sys.exit()
    
    siteNumbers = []
    for i in range(int(S/2)):
        siteNumbers.extend([i+1,i+1])
    
    for i, row in enumerate(SbyS):
        
        name = str()
        if i <= 39: name = 'north'
        else: name = 'south'
        
        if (i+1)%2 > 0:
            SbyS[i] = [0.03, name + str(siteNumbers[i]) + '_all', S] + list(row)
        else:            
            SbyS[i] = [0.03, name + str(siteNumbers[i]) + '_active', S] + list(row)
            
    return SbyS
    
    
    
def getlists(ABV1, ABV2):
    
    """ This function uses mysteries to do magic and return two lists.
    
    ABV1 & ABV2  :  Lists that contain species id#s and the species abundance
                    e.g. ABV1 = [[1, 251], [2, 187], [3, 122], ...]
                    where the first element in each sublist is the species
                    id# and the second element is the abundance.
    """
    
    max1 = ABV1[-1][0]
    max2 = ABV2[-1][0]
    maxID = max([max1, max2])
    
    list1 = [0]*maxID
    list2 = [0]*maxID
    
    for i, sp in enumerate(ABV1):
        
        spID = sp[0]
        ab = sp[1]
        list1[spID-1] = ab 
        
    for i, sp in enumerate(ABV2):
        
        spID = sp[0]
        ab = sp[1]
        list2[spID-1] = ab
        
    return [list1, list2]



def Bray_Curtis(ABV1, ABV2):
    
    """ A function that calculates the Bray-Curtis dissimilarity index, coded
    by Ken Locey.

    Because similarity is more intuitive than dissimilarity, this function
    actually return % similarity.

    A value of 0 means that the two sets being compared have no elements
    in common. A value of 1 means the two sets have identical members.
        
    ABV1 & ABV2  :  Lists that contain species id#s and the species abundance
                    e.g. ABV1 = [[1, 251], [2, 187], [3, 122], ...]
                    where the first element in each sublist is the species
                    id# and the second element is the abundance.
                    
    This function expects the lists to come sorted smallest to largest
    according to species id#  
    """
    list1, list2 = getlists(ABV1, ABV2)
    
    C = float()
    for i, ab in enumerate(list1):
        C += np.abs(ab - list2[i])
    
    C = 100 * (1 - C/(sum(list1) + sum(list2)))
                               
    return C
    


def e_var(SAD):
    
    """ Calculates Smith and Wilson's evenness index. coded by Ken Locey.
    
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
            if i[1] == 'd':
                D += 1
                Ss.append(i[0])
            elif i[1] == 'a':
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
            

    

def reproduction(COM, enVal):
    
    """ A function to simulate reproduction in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    p = randrange(len(COM)) # choose a patch from COM at random    
    if len(COM[p]) == 0: return COM # if the patch is empty, return
    
    i = randrange(len(COM[p])) # choose an individual at random from patch p
    ind = COM[p][i]
    
    if ind[1] == 'd':
        return COM # dormant individuals can't reproduce
    
    optima = 1/ind[0] # the species optima is the reciprocal of its
                      # identity but can also be 1-(1/p), etc.
                
    Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment    
                                 # value of 0 means no disagreement.
                                    
    reproduce = np.random.binomial(1, Err) # higher error means
                                           # low chance of reproducing
    if reproduce == 1:
        COM[p].append(COM[p][i])
        
    return COM
        


def dormancy(COM, enVal):
    
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
    
    p = randrange(len(COM)) # choose a patch from COM at random    
    if len(COM[p]) == 0: return COM # if the patch is empty, return
    
    i = randrange(len(COM[p])) # choose an individual at random from patch p
    ind = COM[p][i]
    
    optima = 1/ind[0] # the species optima is the reciprocal of its
                      # identity but can also be 1-(1/p), etc.
    
    Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment;   
                                 # value of 0 means no disagreement.
                                    
    dormant = np.random.binomial(1, Err) # higher error means higher 
                                         # chance of being dormant
    if dormant == 1:
        ind[1] = 'd' # give the individual a dormant status
    else: 
        ind[1] = 'a' # give the individual an active status
    
    COM[p][i] = ind
    return COM
        



def leaveORdie(COM, enVal):

    """ A function to simulate death/emigration in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    p = randrange(len(COM)) # choose a patch from COM at random    
    if len(COM[p]) == 0: return COM # if the patch is empty, return
    
    i = randrange(len(COM[p])) # choose an individual at random from patch p
    ind = COM[p][i]
    
    optima = 1/ind[0] # the species optima is the reciprocal of its
                      # identity but can also be 1-(1/p), etc.
        
    Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment    
                                 # value of 0 means no disagreement.
                                        
    buhbye = np.random.binomial(1, Err) # higher error means
                                        # higher chance of leaving/death
    if buhbye == 1:
        COM[p].pop(i)
        
    return COM



def disperse(COM1, COM2, enVal):

    """ A function to simulate disperal between two local communities/areas.
    
    COM1 & COM2  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive.    
    """
    
    p1 = randrange(len(COM1))  # choose a patch at random from COM1
    if len(COM1[p1]) == 0: return [COM1, COM2] # if the patch is empty, return
    
    i = randrange(len(COM1[p1])) # choose an individual at random from COM1
    ind = COM1[p1][i]
    
    p2 = randrange(len(COM2))    # choose a patch at random from COM2
    
    optima = 1/ind[0] # the species optima is the reciprocal of its
                      # identity but can also be 1-(1/p), etc.
    
    Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment;    
                                 # value of 0 means no disagreement.
                                        
    dormant = np.random.binomial(1, Err) # higher error means
                                         # higher chance of remaining dormant
    if dormant == 1:
        ind[1] = 'd'   # give the individual a dormant status
    else: 
        ind[1] = 'a'   # give the individual an active status
    
    COM2[p2].append(ind)
    COM1[p1].pop(i)
    
    return COM1, COM2



        
def microbide(imRate, num_patches, lgp, northVal, southVal, time=500):

    """
    imRate  :  number individuals immigrating from the regional pool
               each time step
    
    lgp  :  log-series parameter, typically close to 0.99
    
    state  :  determines whether local communities will have different
              have different environmental values. Value can be:
              'heterogeneous' or 'homogeneous'
            
    time  :  number of time steps to run simulation. Needs to be large
             enough to allow the community to reach a relatively stable state
    """
    
    northCOM = [list([]) for _ in range(num_patches)]
    southCOM = [list([]) for _ in range(num_patches)]
    
    t = 1
    N = imRate
    while t <= time:
        
        """ Immigration from regional pool"""
        propagules = np.random.logseries(lgp, imRate)  # list of propagules
        # the numbers in the list represent species ID numbers
        
        # Note: in the log-series, there are lots of 1's and few very large
        # numbers
        
        for prop in propagules:
            
            northORsouth = np.random.binomial(1, 0.5)
            
            if northORsouth == 0:            
                i = randrange(num_patches)
                northCOM[i].append([prop, 'd']) # add dormant individual
                                                # to north community
            else: 
                i = randrange(num_patches)
                southCOM[i].append([prop, 'd']) # add individual
                                                # to south community
        
        for i in range(N):    
            """ Reproduction """
            northCOM = reproduction(northCOM, northVal)
            southCOM = reproduction(southCOM, southVal)
            
            """ Dormancy """
            northCOM = dormancy(northCOM, northVal)
            southCOM = dormancy(southCOM, southVal)
            
            """ Emigration/Death """
            northCOM = leaveORdie(northCOM, northVal)
            southCOM = leaveORdie(southCOM, southVal)
            
            """ Between patch dispersal """
            northORsouth = np.random.binomial(1, 0.5)
            if northORsouth == 1:    
                northCOM, southCOM = disperse(northCOM, southCOM, southVal)
            else: southCOM, northCOM = disperse(southCOM, northCOM, northVal)
            
        t += 1
        
        N, S, D, A, BP, Evar, nABV = getNSDA(northCOM)
        print 'time:',t, 'north N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        # ABV stands for 'abundance vector'; a list of species identies and
        # their numerical abundances
        
        N, S, D, A, BP, Evar, sABV = getNSDA(southCOM)
        print 'time:',t, 'south N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        BC = Bray_Curtis(nABV, sABV)
        print 'time:',t, 'Bray-Curtis as percent similarity',BC, '\n'
        
    return [northCOM, southCOM]