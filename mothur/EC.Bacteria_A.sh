#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=500gb,walltime=48:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2017_EvolutionCanyon/
module load gcc/4.9.2
module load mothur/1.39.0
mothur EC.Bacteria_A.Batch
qsub EC.Bacteria_B.sh
