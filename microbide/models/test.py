from __future__ import division
import sys, csv, os, multiprocessing
# import EClandscape as land
# import ECMicrobideCore as model
# import ECfunctions as funx

mypath = os.path.dirname(os.path.realpath(__file__))
path = os.path.join(mypath, 'SbyS')
print(path)


cores = multiprocessing.cpu_count()
print('Available Cores = '+str(cores))
