#!/bin/sh

#SBATCH --job-name=J1201

##SBATCH --partition=all-LoPri,all-HiPri # ARGO
#SBATCH --partition=normal,interactive,contrib # HOPPER
#SBATCH --output=/scratch/%u/%x-%j.out # user, jobName, jobID
#SBATCH --error=/scratch/%u/%x-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=%u@gmu.edu

#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=30G

# ARGO
# module load python/3.8.4
# source /home/%u/python-ve/env/bin/activate

module load gnu10/10.3.0-ya
module load python/3.10.1-qb
source /home/sdoan/venv/bin/activate

python ../../../run.py J1201.toml
