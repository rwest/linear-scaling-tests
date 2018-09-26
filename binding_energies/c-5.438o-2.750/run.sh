#!/bin/bash
#SBATCH --job-name=RMGcat
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH -n1
#SBATCH --mem=10Gb
#SBATCH --time=2:00:00

echo $RMGpy
python  $RMGpy/rmg.py -p input.py
