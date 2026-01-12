#!/bin/bash
name=elp_car9_aligned

# note: you need to have water.pdb and elp_car9_aligned.pdb in the same directory where you run the packmol command in line 18 of this script. 

# cp ../structures/water.pdb .

# optional: make the create_packmol.sh executable 
chmod +x ./create_packmol.sh 

python ../../scripts/removeH.py ${name}.pdb

echo '5'|gmx pdb2gmx -f ${name}_noH.pdb -o ${name}_noH.gro -water tip3p -ter

gmx editconf -f ${name}_noH.gro -o ${name}_final.pdb

# this generates the packmol input
./create_packmol.sh 

# run packmol as follows, make sure to change the path to your packmol executable
# this generates a pdb file called elp_car9_aligned_packed.pdb
path-to-packmol/packmol-20.15.1/packmol < packmol_input.inp 

# generates the final .gro file and corresponding topology for downstream simulations

echo '5'| gmx pdb2gmx -f ${name}_packed.pdb -o ${name}_packed_final.gro -water tip3p -ter -ignh

echo '1' 'q'|gmx make_ndx -f ${name}_packed_final.gro -o system.ndx


