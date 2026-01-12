#!/bin/bash
set -euo pipefail
# make an index file that contains a group for protein and water only
echo -e "1|12\nq" | gmx make_ndx -f npt2.gro -o prot_water.ndx

# create a new .mdp file that includes energy groups for protein and water
cp ../../setup_and_run/step3_simulation_protocol/mdp/npt2.mdp ./npt2_prot_water.mdp

echo -e "\n; Energy Groups\nenergygrps              = Protein Water\n" >> npt2_prot_water.mdp

# use the -rerun option to calculate energies for protein and water groups
gmx grompp -f npt2_prot_water.mdp -c npt2.gro -p topol.top -o npt2_prot_water.tpr -n prot_water.ndx -maxwarn 5
gmx mdrun -v -ntmpi 1 -s npt2_prot_water.tpr -rerun npt2.xtc -e npt2_prot_water.edr

# extract protein-protein and protein-water lennard-jones and coulombic energies
echo -e "16 17 20 21"|gmx energy -f npt2_prot_water.edr -o energy.xvg