#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=32,vmem=50gb,walltime=24:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu,kjlocey@indiana.edu
#PBS -m abe
#PBS -j oe
#PBS -o /N/dc2/projects/Lennon_Sequences/2014_EvolutionCanyon/microbide/
module load intel
module load python
cd /N/dc2/projects/Lennon_Sequences/2014_EvolutionCanyon/microbide
python ./models/ECmodel.py >> ECmodel_20140730.log
