export RMGpy=$HOME/Code/RMG-Py
find . -name run.sh -execdir sh -c "sbatch run.sh" \;
