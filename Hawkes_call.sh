#!/bin/bash
#SBATCH --job-name=frog_call
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=250G
#SBATCH --mail-user=eliza.grames@uconn.edu
#SBATCH -o myscript_%j.out
#SBATCH -e myscript_%j.err

echo `hostname`
cd /home/FCAM/egrames/hawkes_vocalizations/
module load R/3.5.1
module load JAGS/4.3.0
Rscript scripts/frog_script.R


