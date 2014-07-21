from __future__ import division
import sys                                            
import numpy as np


""" Core functions for Evolution Canyon Project modeling """

def SpeciesInPatches(COM, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys, SRow2Ys):
    
    nCOM = []
    sCOM = []
    
    NXs = NRowXs.extend(NRowXs)
    NYs = NRow1Ys.extend(NRow2Ys)
    
    SXs = SRowXs.extend(SRowXs)
    SYs = SRow1Ys.extend(SRow2Ys)
    
    NPatches = [list([]) for _ in range(len(NXs))]
    SPatches = [list([]) for _ in range(len(SXs))]
    
    for i, ind in enumerate(COM):
        x, y = ind[2], ind[3]
        
        for p, xcoord in enumerate(NXs):
            if xcoord - .02 <= x and xcoord + .02 >= x:         
                ycoord = NYs[p]
                if ycoord - .02 <= y and ycoord + .02 >= y:
                    NPatches[p].append(ind[0],ind[1])
                    
            elif SXs[p] - .02 <= x and SXs[p] + .02 >= x:         
                ycoord = SYs[p]
                if ycoord - .02 <= y and ycoord + .02 >= y:
                    SPatches[p].append(ind[0],ind[1])
                    
    return [NPatches, SPatches]           
                     
    
    
    
    



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