#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=500gb,walltime=24:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2015_EvolutionCanyon/
module load gcc/4.9.2
module load mothur/1.36.1
mothur EC.Bacteria_A.Batch
qsub EC.Bacteria_B.sh
