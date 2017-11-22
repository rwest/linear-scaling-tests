
# remove (from the repo not the disk) all the solver csv files.
find binding_energies -name solver -execdir sh -c "git rm -r --cached --ignore-unmatch solver" \;
# add just the highest numbered (last) solver csv
find binding_energies -name solver -execdir sh -c "ls solver/*.csv | tail -n 1 | xargs git add" \;
# add everything not ignored in the chemkin directories
find binding_energies -name chemkin -execdir sh -c "git add chemkin" \;
