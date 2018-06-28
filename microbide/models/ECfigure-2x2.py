from __future__ import division
from random import sample
import os
import matplotlib.pyplot as plt
import numpy as np
import sys

mydir = os.path.expanduser('~/GitHub/evolution-canyon/microbide/models')
sys.path.append(mydir)

import ECfunctions as funx


fig = plt.figure()
for i in range(1,5):

    n = 1 * 10**4
    Alpha, Beta = 10, 2

    Nx = np.random.beta(Alpha, Beta, n).tolist()
    #Nx = filter(lambda a: a >= 0.45, Nx)
    Sx = np.random.beta(Beta, Alpha, n).tolist()
    #Sx = filter(lambda a: a <= 0.55, Sx)

    Ny = np.random.uniform(0, 1, len(Nx)).tolist()
    Sy = np.random.uniform(0, 1, len(Sx)).tolist()

    fig.add_subplot(2, 2, i)
    gs = 20

    if i == 1:
        imageS = plt.hexbin(Sx, Sy, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlGn, alpha=0.6)

        imageN = plt.hexbin(Nx, Ny, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlOrBr, alpha=0.6)

    elif i == 2:
        imageN = plt.hexbin(Nx, Ny, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlGn, alpha=0.6)
        imageS = plt.hexbin(Sx, Sy, mincnt=1, gridsize = gs, bins = 'log',
            cmap=plt.cm.YlGn, alpha=0.6)

    if i in [1,2]:
        Ncounts = imageN.get_array()
        Nverts = imageN.get_offsets()
        Nverts = sample(Nverts, 10)
        Nverts = np.array(Nverts)

        Scounts = imageS.get_array()
        Sverts = imageS.get_offsets()
        Sverts = sample(Sverts, 10)
        Sverts = np.array(Sverts)

        NRow1Xs = [0.1] * 10
        NRow2Xs = [0.2] * 10
        NRowYs = [0.1,  0.189,  0.278,  0.367,  0.456, 0.544,  0.633,  0.72,  0.811,  0.9]

        SRow1Xs = [0.8] * 10
        SRow2Xs = [0.9] * 10
        SRowYs = [0.1,  0.189,  0.278,  0.367,  0.456, 0.544,  0.633,  0.72,  0.811,  0.9]

        plt.scatter(NRow1Xs, NRowYs, color='0.9', s=40, edgecolor='0.2')
        plt.scatter(NRow2Xs, NRowYs, color='0.9', s=40, edgecolor='0.2')

        plt.scatter(SRow1Xs, SRowYs, color='0.9', s=40, edgecolor='0.2')
        plt.scatter(SRow2Xs, SRowYs, color='0.9', s=40, edgecolor='0.2')

        plt.xlim(0, 1)
        plt.ylim(0, 1)

        plt.xlabel('Side 1                   Side 2', fontsize=10)

    ''' # optional text
    if i == 1:
        title = 'Birds-eye view of a simulated two-sided landscape with different\n'
        title += 'environments on each side. Each side has 10 plots in 2 rows.'
        plt.text(0.0, 1.1, title, fontsize=10)
    '''

    if i == 3:
        c1 = 'green'
        c2 = 'brown'
    elif i == 4:
        c1 = 'green'
        c2 = 'green'

    if i in [3, 4]:
        DN = funx.get_kdens(Nx)
        plt.plot(DN[0], DN[1], color = c2, lw=3, label= 'Side 1', alpha=0.5)
        plt.fill_between(DN[0], DN[1], 0, color = c2, alpha = 0.3)

        DS = funx.get_kdens(Sx)
        plt.plot(DS[0], DS[1], color = c1, lw=3, label= 'Side 2', alpha=0.5)
        plt.fill_between(DS[0], DS[1], 0, color = c1, alpha = 0.3)

        plt.xlabel('Environmental Condition', fontsize=10)
        plt.ylabel('Frequency', fontsize=10)

        ''' # optional text
        if i == 3:
            txt =  'Below: Distribution of environmental values across sides of the'
            txt += '\nsimulated landscape. The mismatch of species optima to'
            txt += 'local conditions determines\nprobabilities of death,'
            txt += 'reproduction, and transitions to and from dormancy.'

            plt.xlim(0, 1)
            plt.ylim(0, 4.2)
            plt.text(0.0, 4.6, txt, fontsize=10)
        '''

    plt.tick_params(axis='both', labelsize=6)
    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    plt.savefig(mydir+'/figures/final-figs/ECFig-2x2.png',
        dpi=200, bbox_inches='tight', pad_inches=0.1)
