from __future__ import division

import matplotlib # statements for running on IU servers be sure to comment out
matplotlib.use('Agg') # if running on local machine, resets where python looks
from matplotlib import pyplot as plt

import numpy as np
from random import sample

def get_landscape(combo):

    envDiff, enterD, exitD = combo
    n = 1 * 10**3 # size of sample to be draw from a beta distribution

    Alpha, Beta = 10, 3
    Nx = np.random.beta(Alpha, Beta, n).tolist()
    Nx = filter(lambda a: a >= 0.5, Nx)
    Sx = np.random.beta(Beta, Alpha, n).tolist()
    Sx = filter(lambda a: a <= 0.5, Sx)

    Ny = np.random.uniform(0, 1, len(Nx)).tolist()
    Sy = np.random.uniform(0, 1, len(Sx)).tolist()

    if envDiff == 'differ':
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20, cmap=plt.cm.YlOrBr, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20, cmap=plt.cm.YlGn, alpha=0.6)

    elif envDiff == 'same':
        imageN = plt.hexbin(Nx, Ny, mincnt=0, gridsize = 20, cmap=plt.cm.YlGn, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=0, gridsize = 20, cmap=plt.cm.YlGn, alpha=0.6)

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

    return [[NRowYs, NRow1Xs, NRow2Xs, SRowYs, SRow1Xs],
                [SRow2Xs, Ncounts, Nverts, Scounts, Sverts]]
