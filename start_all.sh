export RMGpy=/scratch/westgroup/mazeau/Cat/RMG-Py
find . -name run.sh -execdir sh -c "sbatch run.sh" \;
