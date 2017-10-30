#!/bin/bash
#SBATCH --job-name=RMGcat
#SBATCH --error=error.log
#SBATCH --output=output.log
#SBATCH -n1

echo $RMGpy
python  $RMGpy/rmg.py -p input.py
