#!/bin/sh

#SBATCH --job-name=XID2028_nirspec_ifu

#SBATCH --partition=all-LoPri,all-HiPri # ARGO
##SBATCH --partition=normal # HOPPER
#SBATCH --output=/scratch/%u/%x-%j.out # user, jobName, jobID
#SBATCH --error=/scratch/%u/%x-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=sdoan2@gmu.edu

#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH --mem=30G


#module load python/3.8.6-ff # HOPPER
module load python/3.8.4 # ARGO
source /home/sdoan2/python-ve/env/bin/activate

python ../../run_target.py XID2028_nirspec_ifu.toml
