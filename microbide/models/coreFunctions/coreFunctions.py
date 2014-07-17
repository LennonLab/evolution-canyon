from __future__ import division
import sys                                            
import numpy as np
import random
from random import randrange, choice
import decimal


""" Core functions for Evolution Canyon Project modeling """

def get_param(ID, param_dict):
    """ A function to set species specific parameter values """
    param = np.random.uniform(0.01, 0.9)
    param_dict[ID] = param
    
    return param_dict
    

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
            

    

def reproduction(COM):
    
    """ A function to simulate reproduction in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    for i, patch in enumerate(COM):
        while patch:    
            ind = randrange(len(patch)) # choose an individual at random from patch p
            status = patch[ind][1]
            
            if status == 0:
                COM[i].append(patch[ind])
                break
                                  
    return COM
        


def dormancy(COM, enVal, optima_dict):
    
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
    i = randrange(len(COM[p])) # choose an individual at random from patch p
    
    ind = COM[p][i]
    optima = optima_dict[ind[0]] # the species environmental optima 
    
    Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment;   
                                 # value of 0 means no disagreement.
                                    
    x = np.random.binomial(1, Err) # higher error means higher 
                                   # chance of being dormant
    
    if ind[1] == 0: COM[p][i][1] = x
    elif x == 0: COM[p][i][1] = 0
    elif x == 1: COM[p][i][1] += 1
    
    return COM
        



def death_emigration(COM, enVal, optima_dict, ceiling=50):

    """ A function to simulate death/emigration in a local community.
    
    COM  :  the local community; a vector of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    enVal  :  the environmental variable value of the local community. ranges 
              between 0.0 and 1.0, inclusive. """
    
    for index, patch in enumerate(COM):
        while len(patch) > ceiling:    
            
            i = randrange(len(patch)) # choose an individual at random from patch p
            ind = patch[i]
    
            optima = optima_dict[ind[0]] # the species environmental optima         
            Err = np.abs(optima - enVal) # difference between species
                                 # optima and patch environment    
                                 # value of 0 means no disagreement.
                                        
            buhbye = np.random.binomial(1, Err) # higher error means
                                        # higher chance of leaving/death
            if buhbye == 1:
                patch.pop(i)
        
        COM[index] = patch
        
    return COM



def disperse(COM1, COM2, enVal, between=True):

    """ A function to simulate disperal between two local communities/areas.
    
    COM1 & COM2  :  the local communities; vectors of lists where each list is a patch
            and each patch contains smaller lists that represent the information
            for individuals, e.g. COM[0][0] is the zeroth individual of the 
            zeroth patch of COM.
    
    """
    if between == True:
        p = randrange(len(COM1))  # choose a patch at random from COM1
        i = randrange(len(COM1[p])) # choose an individual at random from COM1
        ind = COM1[p].pop(i)
    
        p = randrange(len(COM2)) # choose a patch at random from COM2
        COM2[p].append(ind)

    elif between == False:
        p = randrange(len(COM1))  # choose a patch at random from COM1
        i = randrange(len(COM1[p])) # choose an individual at random from COM1
        ind = COM1[p].pop(i)
    
        p = randrange(len(COM1)) # choose a patch at random from COM2
        COM1[p].append(ind)
        
    return COM1, COM2




def immigration(COMs, optima_dict, env_optima, decay_dict, dorm_decay, im):
    
    for i, COM in enumerate(COMs):
        
        for j, site in enumerate(COM):     
            propagules = np.random.logseries(0.9, im)  # list of propagules 
    
            for prop in propagules:
            
                if prop not in optima_dict:
                    optima_dict = get_param(prop, optima_dict)
            
                if prop not in decay_dict:
                    decay_dict = get_param(prop, decay_dict)
                
                COMs[i][j].append([prop, choice([0, 1])]) # 0 is active, 1 is dormant
                                               
                                                
    return COMs
              
                            
    
        
def microbide(num_patches, lgp, northVal, southVal, env_optima, dorm_decay, btn):

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
    time=200
    ceil = 50 #     number of individuals per patch
    
    optima_dict = {}
    decay_dict = {}
    northCOM = [list([]) for _ in range(num_patches)]
    southCOM = [list([]) for _ in range(num_patches)]
    northCOM, southCOM = immigration([northCOM, southCOM], optima_dict,
                                    env_optima, decay_dict, dorm_decay, ceil)
    
    
    im = 1 #     number of individuals immigrating per unit time
    t = 1
    while t <= time:
        
        """ Immigration from regional pool. At time = 1, the community patches
        are populated with the same large number of individuals. Then a zero-sum
        dynamic is enforced where death, immigration, emigration, etc. occur but
        never change the abundance in a patch. The zero-sum dynamic is our way
        of keeping patch abundance from exploding and/or drastically shrinking.
        
        Species IDs of immigrants are chosen at random from a log-series
        distribution, i.e., a monotonically decreasing frequency distribution.
        
        """
        
        """ Immigration """
        northCOM, southCOM = immigration([northCOM, southCOM], optima_dict,
                                    env_optima, decay_dict, dorm_decay, im)
    
        """ Reproduction """
        northCOM = reproduction(northCOM)
        southCOM = reproduction(southCOM)
            
        """ Dormancy """
        northCOM = dormancy(northCOM, northVal, optima_dict)
        southCOM = dormancy(southCOM, southVal, optima_dict)
        
        """ Between patch dispersal """
        northCOM, southCOM = disperse(northCOM, southCOM, southVal, btn)
        southCOM, northCOM = disperse(southCOM, northCOM, northVal, btn)    
            
        """ Emigration & Death """
        northCOM = death_emigration(northCOM, northVal, optima_dict, ceil)
        southCOM = death_emigration(southCOM, southVal, optima_dict, ceil)
            
        t += 1
        
        N, S, D, A, BP, Evar, nABV = getNSDA(northCOM)
        print 'time:',t, 'north N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        # ABV stands for 'abundance vector'; a list of species identies and
        # their numerical abundances
        
        N, S, D, A, BP, Evar, sABV = getNSDA(southCOM)
        print 'time:',t, 'south N:',N,'S:',S,'D:',D,'A:',A,'BP:',BP,'Evar:',Evar
        
        BC = Bray_Curtis(nABV, sABV)
        print 'time:',t, 'Bray-Curtis as percent similarity',BC, '\n'
    
    #print 'number of sites in north and south:',len(northCOM), len(southCOM)    
    return [northCOM, southCOM]