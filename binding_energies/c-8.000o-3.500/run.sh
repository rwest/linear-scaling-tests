#!/bin/bash
#SBATCH --job-name=RMGcat
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH -n1
#SBATCH --mem=120Gb
#SBATCH --time=24:00:00

echo $RMGpy
python  $RMGpy/rmg.py -p input.py
