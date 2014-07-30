from __future__ import division 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np
from random import sample

def get_landscape(combo):
    
    envDiff, enterD, exitD = combo                
    n = 1 * 10**4 # size of sample to be draw from a beta distribution
    
    #envDiff = 'differ' # These two statement are used for informal testing
    #envDiff = 'same'
    
    Alpha, Beta = 10, 2
    Nx = np.random.beta(Alpha, Beta, n).tolist()
    Nx = filter(lambda a: a >= 0.5, Nx)
    Sx = np.random.beta(Beta, Alpha, n).tolist()
    Sx = filter(lambda a: a <= 0.5, Sx)
        
    Ny = np.random.uniform(0, 1, len(Nx)).tolist()
    Sy = np.random.uniform(0, 1, len(Sx)).tolist()
    
    if envDiff == 'differ':
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20,
            bins = 'log', cmap=plt.cm.YlOrBr, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20,
            bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
    
    elif envDiff == 'same':
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20,
            bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20,
            bins = 'log', cmap=plt.cm.YlGn, alpha=0.6)
                
    Ncounts = imageN.get_array()
    Nverts = imageN.get_offsets()
    Nverts = sample(Nverts, 10)
    Nverts = np.array(Nverts)
    
    Scounts = imageS.get_array()
    Sverts = imageS.get_offsets()
    Sverts = sample(Sverts, 10)
    Sverts = np.array(Sverts)
    
    NRow1Xs = [0.06] * 10  
    NRow2Xs = [0.12] * 10
    NRowYs = [0.26, 0.32, 0.38, 0.44, 0.50, 0.56, 0.62, 0.68, 0.74, 0.80] 
    								
    SRow1Xs = [0.88] * 10  
    SRow2Xs = [0.94] * 10
    SRowYs = [0.26, 0.32, 0.38, 0.44, 0.50, 0.56, 0.62, 0.68, 0.74, 0.80] 
    								
    return [[NRowYs, NRow1Xs, NRow2Xs, SRowYs, SRow1Xs],
                [SRow2Xs, Ncounts, Nverts, Scounts, Sverts]]
                
