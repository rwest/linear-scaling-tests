export RMGpy=/gss_gpfs_scratch/westgroup/Cat/RMG-Py
find . -name run.sh -execdir sh -c "sbatch run.sh" \;
