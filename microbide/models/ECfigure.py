from __future__ import division 
from random import sample

#import matplotlib  # These two statements are for running in IU servers
#matplotlib.use('Agg')

import ECfunctions as funx
import matplotlib.pyplot as plt 
import numpy as np
import sys

path = '/Users/lisalocey/Desktop/evolution-canyon/microbide/models/'

def ECfig(envDiff):
    
    n = 1 * 10**4
    Alpha, Beta = 10, 3
    
    Nx = np.random.beta(Alpha, Beta, n).tolist()
    #Nx = filter(lambda a: a >= 0.5, Nx)
    Sx = np.random.beta(Beta, Alpha, n).tolist()
    #Sx = filter(lambda a: a <= 0.5, Sx)
        
    Ny = np.random.uniform(0, 1, len(Nx)).tolist()
    Sy = np.random.uniform(0, 1, len(Sx)).tolist()
    
    fig = plt.figure()
    fig.add_subplot(1, 1, 1)
    gs = 20
    
    if envDiff == 'differ':
        # extra hexbin layer to beautify the landscape plot
        #plt.hexbin(Nx, Ny, mincnt=0, cmap=plt.cm.jet, gridsize = gs, alpha=0.4)
        #plt.hexbin(Sx, Sy, mincnt=0, cmap=plt.cm.jet, gridsize = gs, alpha=0.4)
    
        imageN = plt.hexbin(Nx, Ny, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlOrBr, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlGn, alpha=0.6)
    
    elif envDiff == 'same':
        # extra hexbin layer to beautify the landscape plot
        #plt.hexbin(Nx, Ny, mincnt=0, cmap=plt.cm.jet, gridsize = gs, alpha=0.4)
        #plt.hexbin(Sx, Sy, mincnt=0, cmap=plt.cm.jet, gridsize = gs, alpha=0.4)
        
        imageN = plt.hexbin(Nx, Ny, mincnt=1, gridsize = gs, bins = 'log', 
            cmap=plt.cm.YlGn, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=1, gridsize = gs, bins = 'log', 
            cmap=plt.cm.YlGn, alpha=0.6)
                
    Ncounts = imageN.get_array()
    Nverts = imageN.get_offsets()
    Nverts = sample(Nverts, 10)
    Nverts = np.array(Nverts)
    
    Scounts = imageS.get_array()
    Sverts = imageS.get_offsets()
    Sverts = sample(Sverts, 10)
    Sverts = np.array(Sverts)
    
    NRow1Xs = [0.15] * 10  
    NRow2Xs = [0.25] * 10
    NRowYs = [0.26, 0.32, 0.38, 0.44, 0.50, 0.56, 0.62, 0.68, 0.74, 0.80] 
    
    SRow1Xs = [0.75] * 10  
    SRow2Xs = [0.85] * 10
    SRowYs = [0.26, 0.32, 0.38, 0.44, 0.50, 0.56, 0.62, 0.68, 0.74, 0.80]
    
    plt.scatter(NRow1Xs, NRowYs, color='w', s=40, alpha=0.8)							
    plt.scatter(NRow2Xs, NRowYs, color='w', s=40, alpha=0.8)							
    
    plt.scatter(SRow1Xs, SRowYs, color='w', s=40, alpha=0.8)							
    plt.scatter(SRow2Xs, SRowYs, color='w', s=40, alpha=0.8)							
    
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    
    plt.ylabel('West to East', fontsize=14)
    plt.xlabel('South                                               North', fontsize=14)
    
    if envDiff == 'differ':
        title = 'Birds-eye view of the environmental landscape related to mesic (green)'
        title += '\n& xeric (brown) conditions. Each side has 10 plots in 2 rows.'
        plt.title(title, fontsize=13)

    if envDiff == 'same':
        title = 'Birds-eye view of the environmental landscape related to mesic (green)'
        title += '\n conditions. Each side of the simulated area has 10 plots in 2 rows.'
        plt.title(title, fontsize=13)
    
    
    plt.savefig(path+'EC_Heat'+envDiff+'.png', 
        dpi=600,bbox_inches='tight',pad_inches=0.1)
    
    # A SECOND FIGURE
    fig = plt.figure()
    fig.add_subplot(1, 1, 1)
    
    if envDiff == 'differ':
        DN = funx.get_kdens(Nx)
        plt.plot(DN[0], DN[1], color = 'brown', lw=3, label= 'North', alpha=0.5)
        plt.fill_between(DN[0], DN[1], 0, color = 'brown', alpha = 0.3)
                
    DS = funx.get_kdens(Sx)
    plt.plot(DS[0], DS[1], color = 'green', lw=3, label = 'South', alpha=0.5)
    plt.fill_between(DS[0], DS[1], 0, color = 'green', alpha = 0.3)
    
    plt.xlabel('Environmental Condition', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    
    txt =  'Distribution(s) of environmental values across patches. '
    txt += 'Mismatch\nof species optima influences the probability of '
    txt += 'reproduction,\ndeath, and transitions to and from dormancy.'
    
    plt.xlim(0, 1)
    plt.ylim(0, 8)
    plt.text(0.05, 5.2, txt, fontsize=12)
    
    plt.savefig(path+'figures/ECkdens_'+envDiff+'.png', 
        dpi=600,bbox_inches='tight',pad_inches=0.1)
    
    plt.subplots_adjust(wspace=0.0, hspace=0.4)
    plt.show()    
    
    return #[NRowYs, NRow1Xs, NRow2Xs, SRowYs, SRow1Xs,
            #    SRow2Xs, Ncounts, Nverts, Scounts, Sverts]
                

ECfig(envDiff='differ')