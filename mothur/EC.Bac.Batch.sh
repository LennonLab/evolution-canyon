#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=32,vmem=500gb,walltime=96:00:00
#PBS -M lennonj@indiana.edu 
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/EvolutionCanyon3
module load gcc
module load mothur/1.31.2 
mothur EC.Bacteria.Batch
