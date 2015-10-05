from __future__ import division
import sys, csv, os, time
import multiprocessing # Python tool for parallelizing
import EClandscape as land
import ECMicrobideCore as model
import ECfunctions as funx

mypath = os.path.dirname(os.path.realpath(__file__))
path = os.path.join(os.path.split(mypath)[0], 'SbyS')
cores = multiprocessing.cpu_count()

print('Output Path = '+str(path))
print('Available Cores = '+str(cores))
