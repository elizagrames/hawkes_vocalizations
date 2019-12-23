#!/bin/bash
#SBATCH --job-name=Hawkes_process
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=30G
#SBATCH --mail-user=eliza.grames@uconn.edu
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err

echo `hostname`
module load R/3.5.1
module load JAGS/4.3.0
cd hawkes_vocalizations/scripts/
Rscript Hawkes_script_for_paper.R

