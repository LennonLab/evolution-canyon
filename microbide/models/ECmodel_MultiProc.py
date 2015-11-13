from __future__ import division
import sys, csv, os, time, datetime
import multiprocessing # Python tool for parallelizing
import EClandscape as land
import ECMicrobideCore as model
import ECfunctions as funx

mypath = os.path.dirname(os.path.realpath(__file__))
path = os.path.join(os.path.split(mypath)[0], 'SbyS')
cores = multiprocessing.cpu_count()

print(datetime.date.today())
print('Output Path = '+str(path))
print('Available Cores = '+str(cores))

###########################  GET CONDITIONS  ###################################

""" Code to run the microbide model and generates site-by-species matrices."""

num_patches = 20 # number of patches on each side of Evolution Canyon (EC)
lgp = 0.99 # log-series parameter; underlying structure of regional pool
N = 2 * 10**7 # Starting total abundance across the landscape
T = 10**6 # Time parameter

conditions = [[0, 'same', 'rand', 'rand'],
             [1, 'differ', 'rand', 'rand'],
             [2, 'same',  'env',  'env'],
             [3, 'differ', 'rand', 'env'],
             [4, 'differ', 'env',  'rand'],
             [5, 'differ',  'env', 'env']]

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

####################  GENERATE SITE BY SPECIES DATA  ###########################
def worker(combo):

    condition, envDiff, enterD, exitD = combo

    combo.pop(0) # removes first element of combo (i.e., condition)

    landscapeLists = land.get_landscape(combo) # characterizing the landscape

    NRowXs, NRow1Ys, NRow2Ys, SRowXs, SRow1Ys = landscapeLists[0]
    SRow2Ys, Ncounts, Nverts, Scounts, Sverts = landscapeLists[1]

    COM = model.microbide(combo, Ncounts, Nverts, Scounts, Sverts, N, T, condition)
                            # run the model & return the communities

    nCOM, sCOM = funx.SpeciesInPatches(COM, NRowXs, NRow1Ys, NRow2Ys,
                                            SRowXs, SRow1Ys, SRow2Ys)


    SbyS = funx.get_SitebySpecies([nCOM, sCOM]) # get site by species matrix
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

    fileName = os.path.join(path, 'Condition'+str(condition+1))
    OUT = open(fileName + '.txt','w')

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

    return

if __name__ == '__main__':
    pool = multiprocessing.Pool()

    # Parallel map
    tic = time.time()
    results = pool.map(worker, conditions)
    toc = time.time()

    # pool.join()
    # pool.close()

    # Serial map
    #tic2 = time.time()
    #results = map(worker, conditions)
    #toc2 = time.time()

    print 'Parallel processing time: ', str(round((toc-tic)/60, 1)), ' minutes'
