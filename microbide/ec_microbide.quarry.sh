#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,walltime=10:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu,kjlocey@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2014_EvolutionCanyon/microbide
python ./models/ECmodel.py >> ECmodel_N200K_t2M.log
