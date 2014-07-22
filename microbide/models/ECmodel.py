from __future__ import division 
from random import choice, randrange, sample
import matplotlib.pyplot as plt 
import numpy as np
import scipy as sc
from scipy import stats
from scipy.stats import gaussian_kde, sem
import sys, csv

#sys.path.append("/Users/lisalocey/Desktop/evolution-canyon/microbide/models/coreFunctions")
#import coreFunctions as cf


def get_kdens(summands): # Finds the kernel density function across a sample
    
    density = gaussian_kde(summands)
    #xs = np.linspace(float(min(summands)),float(max(summands)), len(summands))
    xs = np.linspace(0,1, len(summands))
    density.covariance_factor = lambda : 0.1
    density._compute_covariance()
    
    return [xs, density(xs)]
    
    
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
    
    
def get_match(Nverts, Sverts, x1, y1, Ncounts, Nmax, Scounts, Smax, opt): 

    dmin = 10**6
    Nenv = 0
        
    for j, vert in enumerate(Nverts):
        x2, y2 = vert
            
        dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if dist < dmin:
            dmin = dist
            Nenv = np.abs(opt - Ncounts[j]/Nmax) # mismatch to North niche
    
    dmin = 10**6
    Senv = 0
    
    for j, vert in enumerate(Sverts):
        x2, y2 = vert
        
        dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if dist < dmin:
            dmin = dist
            Senv = np.abs(opt - Scounts[j]/Smax) # mismatch to South niche
    
    match = 1 - (Senv * Nenv)
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
        


def SpeciesInPatches(COM, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys, SRow2Ys):
    
    NXs = NRowXs + NRowXs
    NYs = NRow1Ys + NRow2Ys
    
    SXs = SRowXs + SRowXs
    SYs = SRow1Ys + SRow2Ys
    
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
      
    
    
def ECfig(combo):
    
    envDiff, enterD, exitD = combo    
                
    n = 1 * 10**4
    envDiff = 'differ'
    envDiff = 'same'
    
    Alpha, Beta = 10, 2
    Nx = np.random.beta(Alpha, Beta, n).tolist()
    Nx = filter(lambda a: a >= 0.5, Nx)
    Sx = np.random.beta(Beta, Alpha, n).tolist()
    Sx = filter(lambda a: a <= 0.5, Sx)
        
    Ny = np.random.uniform(0, 1, len(Nx)).tolist()
    Sy = np.random.uniform(0, 1, len(Sx)).tolist()
    
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    
    if envDiff == 'differ':
        plt.hexbin(Nx, Ny, mincnt=0, cmap=plt.cm.jet, gridsize = 20, alpha=0.4)
        plt.hexbin(Sx, Sy, mincnt=0, cmap=plt.cm.jet, gridsize = 20, alpha=0.4)
    
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20, bins = 'log', cmap=plt.cm.YlOrBr, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20, bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
    
    elif envDiff == 'same':
        plt.hexbin(Nx, Ny, mincnt=0, cmap=plt.cm.jet, gridsize = 20, alpha=0.4)
        plt.hexbin(Sx, Sy, mincnt=0, cmap=plt.cm.jet, gridsize = 20, alpha=0.4)
        
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20, bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20, bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
                
    Ncounts = imageN.get_array()
    Nverts = imageN.get_offsets()
    Nverts = sample(Nverts, 10)
    Nverts = np.array(Nverts)
    
    Scounts = imageS.get_array()
    Sverts = imageS.get_offsets()
    Sverts = sample(Sverts, 10)
    Sverts = np.array(Sverts)
    
    NRowXs = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.4] 
    NRow1Ys = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55] 
    NRow2Ys = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85] 
    								
    for i, val in enumerate(NRowXs):
        plt.plot(NRowXs[i], NRow1Ys[i],'ms', markersize=15, zorder=100, alpha=0.5)
        plt.plot(NRowXs[i], NRow2Ys[i],'ms', markersize=15, zorder=100, alpha=0.5)
  		
    SRowXs = [0.62, 0.66, 0.7, 0.74, 0.78, 0.82, 0.86, 0.90, 0.94, 0.98] 
    SRow1Ys = list(NRow1Ys)
    SRow1Ys.reverse()
    SRow2Ys = list(NRow2Ys)
    SRow2Ys.reverse()
    			
    					
    for i, val in enumerate(SRowXs):
        plt.plot(SRowXs[i], SRow1Ys[i],'ms', markersize=15, zorder=100, alpha=0.5)
        plt.plot(SRowXs[i], SRow2Ys[i],'ms', markersize=15, zorder=100, alpha=0.5)
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    
    plt.ylabel('West to East', fontsize=16)
    plt.xlabel('South                                               North', fontsize=14)
    
    if envDiff == 'differ':
        title = 'Birds-eye view of the environmental landscape related to mesic (green)'
        title += '\n& xeric (brown) conditions. Each side has 10 plots in 2 rows.'
        plt.title(title, fontsize=13)

    if envDiff == 'same':
        title = 'Birds-eye view of the environmental landscape related to mesic (green)'
        title += '\n conditions. Each side of the simulated area has 10 plots in 2 rows.'
        plt.title(title, fontsize=13)
    
            
    # A SECOND PLOT
    ax = fig.add_subplot(2, 1, 2)
    
    if envDiff == 'differ':
        DN = get_kdens(Nx)
        plt.plot(DN[0], DN[1], color = 'brown', lw=3, label= 'North', alpha=0.5)
        plt.fill_between(DN[0], DN[1], 0, color = 'brown', alpha = 0.3)
                
    DS = get_kdens(Sx)
    plt.plot(DS[0], DS[1], color = 'green', lw=3, label = 'South', alpha=0.5)
    plt.fill_between(DS[0], DS[1], 0, color = 'green', alpha = 0.3)
    
    plt.xlabel('Environmental optima', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    
    txt =  'Distribution(s) of environmental optima among species. A species\'\n'
    txt += 'optima determines the individual probaiblity of reproduction,\n'
    txt += 'death, and transitions into and out of dormancy.'
    
    plt.xlim(0, 1)
    plt.ylim(0, 8)
    plt.text(0.05, 5.2, txt, fontsize=12)
    
    #plt.savefig('/Users/lisalocey/Desktop/EC_'+envDiff+'.png', 
    #dpi=600,bbox_inches='tight',pad_inches=0.1)
    
    #plt.subplots_adjust( wspace=0.0, hspace=0.4)
    #plt.show()    
    
    #sys.exit()
    return [fig, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys,
                SRow2Ys, Ncounts, Nverts, Scounts, Sverts]
    
fig = plt.figure()


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


######################### COMMUNITY SIMULATION FUNCTION ########################
def microbide(combo, fig, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys,
                            SRow2Ys, Ncounts, Nverts, Scounts, Sverts, ic):    
    
    envDiff, enterD, exitD = combo 
    
    Nmax = max(Ncounts) # max heat value
    Smax = max(Scounts) # max heat value
    
    oDict, dDict, COM = [ {}, {}, [] ]
    COM = immigration(COM, oDict, dDict, 2*10**6)
    
    for t in range(10**5):
        
        """ Immigration """
        x = np.random.binomial(1, 0.2)
        if x == 1:
            COM = immigration(COM, oDict, dDict, 1)
    
    
        """ Death """
        i = randrange(len(COM))
        spID, state, x1, y1, opt = [COM[i][0], COM[i][1], COM[i][2], 
                                            COM[i][3], oDict[COM[i][0]]]
        
        match = get_match(Nverts, Sverts, x1, y1, Ncounts,
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
        
        match = get_match(Nverts, Sverts, x1, y1, Ncounts,
                                    Nmax, Scounts, Smax, opt)
        
        x = np.random.binomial(1, match)
        if state==0 and x==1: COM.append(COM[i]) # reproduce from a high match   
	
	    
	""" Transition between activity and dormancy """
	i = randrange(len(COM))
        ind = COM[i]
        spID, state, x1, y1, opt = [ind[0], ind[1], ind[2], 
                                            ind[3], oDict[ind[0]]]
        
        if state == 0: # if active
            match = get_match(Nverts, Sverts, x1, y1, Ncounts, 
                                    Nmax, Scounts, Smax, opt)
            
            if enterD == 'random':
                match1 = 0
	        while match1 <=0 or match1 >=1: 
	           match1 = np.random.normal(match, 0.1, 1)
	        match = float(match1)
	    
	    x = np.random.binomial(1, match)
	    if x == 0: COM[i][1] = 1 # go dormant with a low match              
	
	if state == 1: # if dormant
	    match = get_match(Nverts, Sverts, x1, y1, Ncounts, 
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

        
    
####################  GENERATE SITE BY SPECIES DATA  ###########################
for ic, combo in enumerate(conditions):
    
    envDiff, enterD, exitD = combo 
    fig, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys, SRow2Ys, Ncounts, Nverts, Scounts, Sverts =  ECfig(combo) # characterizing landscape

    COM = microbide(combo, fig, NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys,
                            SRow2Ys, Ncounts, Nverts, Scounts, Sverts, ic)
                            # run the model & return the communities
    
    nCOM, sCOM = SpeciesInPatches(COM, NRowXs, NRow1Ys, NRow2Ys,
                                            SRowXs, SRow1Ys, SRow2Ys)                       
    
    SbyS = get_SitebySpecies([nCOM, sCOM]) # get site by species matrix 
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
                                        
    path = '/Users/lisalocey/Desktop/evolution-canyon/microbide/SbyS/'
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