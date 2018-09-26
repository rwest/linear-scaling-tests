#!/bin/bash
#SBATCH --job-name=RMGcat
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH -n1
#SBATCH --mem=25Gb
#SBATCH --time=5:00:00

echo $RMGpy
python  $RMGpy/rmg.py -p input.py
