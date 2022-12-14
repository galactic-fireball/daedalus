#!/bin/sh

#SBATCH --job-name=IIZw096_RATIO

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

echo "$LINE"
#python ../run_program.py --config IIZw096_pipeline_1.8.4.toml --download all
#python ../run_program.py --config IIZw096_pipeline_1.8.4.toml --pipeline
#python ../run_program.py --config IIZw096_pipeline_1.8.4.toml --cube
#python ../run_program.py --config IIZw096_pipeline_1.8.4.toml --line $LINE
python ../run_program.py --config IIZw096_pipeline_1.8.4.toml --ratio
