#!/bin/bash
#SBATCH --job-name=slurm_sample_
#SBATCH --mail-type=ALL
#SBATCH --mail-user=andrew.gibbs@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 
#SBATCH --mem=4gb
#SBATCH --time=120:00:00
#SBATCH --array=1-10
#SBATCH --output=slurm_sample__%A-%a.log

pwd; hostname; date
module load python

python /home/andrew.gibbs/paul.torrey/andrew.gibbs/Research/CGM_Fitting/slurm_array.py $SLURM_ARRAY_TASK_ID 

date



