from __future__ import division
from random import choice, randrange
import numpy as np
from random import shuffle
import sys
import ECfunctions as funx


def immigration(envDiff, COM, oDict1, oDict2, dDict, im, lgp):

    propagules = np.random.logseries(lgp, im) # list of propagules (was 0.99, 0.92 MEM)
    for p in propagules:
        if p not in oDict1:

            oDict1[p] = np.random.uniform(0.0001, 0.9999, 1).tolist()[0]
            if envDiff == 'same':
                oDict2[p] = float(oDict1[p])

            elif envDiff == 'differ':
                oDict2[p] = np.random.uniform(0.0001, 0.9999, 1).tolist()[0]

            dDict[p] =  np.random.uniform(0, 0.01, 1).tolist()[0]

        x = np.random.uniform(0.01, 0.99, 1).tolist()[0]
        y = np.random.uniform(0.01, 0.99, 1).tolist()[0]

        COM.append([p, choice([0, 1]), x, y]) # 0 is active, 1 is dormant

    return COM


######################### COMMUNITY SIMULATION FUNCTION ########################
def microbide(combo, Ncounts, Nverts, Scounts, Sverts, N, T, ic, lgp, im, dkern):

    envDiff, enterD, exitD = combo

    Nmax = max(Ncounts) # max heat value
    Smax = max(Scounts) # max heat value

    oDict1, oDict2, dDict, COM = [ {}, {}, {}, [] ]

    # N = 2 * 10**6 # Starting total abundance across the landscape
    COM = immigration(envDiff, COM, oDict1, oDict2, dDict, N, lgp)

    for t in range(T):

        procs = range(5)
        shuffle(procs)

        for proc in procs:
            if proc == 0 or t == 0:

                """ Immigration """
                # im = probability of an individual immigrating, 0.0 means that

                x = np.random.binomial(1, im)
                if x == 1:
                    COM = immigration(envDiff, COM, oDict1, oDict2, dDict, x, lgp)

            elif proc == 1:
                """ Death """
                i = randrange(len(COM))
                ind = COM[i]

                spID, state, x1, y1 = [ind[0], ind[1], ind[2], ind[3]]
                opt1, opt2 =  oDict1[ind[0]], oDict2[ind[0]]

                match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                                    Nmax, Scounts, Smax, opt1, opt2)

                if state == 0: # if the individual is active
                    x = np.random.binomial(1, match)
                    if x == 0: COM.pop(i) # likely to die from a low match

       	        elif state == 1: # if the individual is dormant
                    d = dDict[spID] # the individual's dormancy parameter serves as a
                            # general "strength of dormant response" variable

                    x = np.random.binomial(1, d * match) # multiplicative interaction
                                    # between dormancy decay and mismatch to env.
                    if x == 0:
                        x = np.random.binomial(1, d)
                        if x == 1:
                            COM.pop(i) # death in dormancy, inherently becomes a
                            # negative binomial distributed random variable, i.e.,
                            # the probability of dying in dormancy after x resamples

            elif proc == 2:
                """ Reproduction """
                i = randrange(len(COM))
                ind = COM[i]

                spID, state, x1, y1 = [ind[0], ind[1], ind[2], ind[3]]
                opt1, opt2 =  oDict1[ind[0]], oDict2[ind[0]]

                match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                                    Nmax, Scounts, Smax, opt1, opt2)

                x = np.random.binomial(1, match)
                if state == 0 and x == 1:
                    COM.append(COM[i]) # reproduce from a high match
	           #print 'reproduce with match of,', match

            elif proc == 3:
                """ Transition between activity and dormancy """

                i = randrange(len(COM))
                ind = COM[i]

                spID, state, x1, y1 = [ind[0], ind[1], ind[2], ind[3]]
                opt1, opt2 =  oDict1[ind[0]], oDict2[ind[0]]

                if state == 0: # if active
                    if enterD == 'env':
       	                match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                            Nmax, Scounts, Smax, opt1, opt2)

                    elif enterD == 'rand':
                        match = 0.5 # going super simple. everybody gets 0.5

                    x = np.random.binomial(1, match)
                    if x == 0:
	               COM[i][1] = 1 # go dormant with a low match
	               #print 'go dormant with a match of', round(match,2)

                elif state == 1: # if dormant

                    if exitD == 'env':
       	                match = funx.get_match(Nverts, Sverts, x1, y1, Ncounts,
                            Nmax, Scounts, Smax, opt1, opt2)

                    elif exitD == 'rand':
	               match = 0.5  # going super simple. everybody gets 0.5

                x = np.random.binomial(1, match)
                if x == 1:
                    COM[i][1] = 0 # activate with a high match
                    #print 'go active with a match of', round(match,2)

   	    elif proc == 4:
                """ Dispersal """
                i = randrange(len(COM))
                x = np.random.normal(loc=COM[i][2], scale=dkern)
                y = np.random.normal(loc=COM[i][3], scale=dkern)

                if x > 1 or x < 0: COM.pop(i)
                elif y > 1 or y < 1: COM.pop(i)

        t += 1
        if t%1000 == 0:
            Slist = []
            for i, ind in enumerate(COM):
                Slist.append(ind[0])

            print 'condition',ic+1,' lgp:',round(lgp,3),' time:',t, 'N:', len(COM), "S:", len(list(set(Slist)))

    return COM
