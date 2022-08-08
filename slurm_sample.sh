#!/bin/bash
#SBATCH --job-name=slurm_sample_
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.gibbs@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=4gb
#SBATCH --time=120:00:00
#SBATCH --output=slurm_sample_.log

pwd; hostname; date
module load python

python /home/andrew.gibbs/paul.torrey/andrew.gibbs/Research/CGM_Fitting/slurm_sample.py

date



