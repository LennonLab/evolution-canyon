#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
#PBS -o /N/dc2/projects/Lennon_Sequences/2016_EvolutionCanyon/
module load gcc/4.9.2
module load intel
module load python
now=$(date +"%Y%m%d")
cd /N/dc2/projects/Lennon_Sequences/2016_EvolutionCanyon/
python ./microbide/ECmodel.py >> ./ECmodel_$now.log
