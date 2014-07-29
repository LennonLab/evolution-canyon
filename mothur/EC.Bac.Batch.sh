#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=32,vmem=500gb,walltime=96:00:00
#PBS -M mmuscare@indiana.edu,lennonj@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/2014_EvolutionCanyon/EvolutionCanyon_Silva
module load gcc
module load mothur/1.32.1 
mothur EC.Bacteria.Batch
