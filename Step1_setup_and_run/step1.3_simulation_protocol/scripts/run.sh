#!/bin/bash
name=elp_car9_aligned
mdp_dir=../mdp
# counterions
gmx grompp -f ${mdp_dir}/ions.mdp -c ${name}_packed_final.gro -p topol.top -n system.ndx -o ions.tpr -maxwarn 3

echo "SOL"| gmx genion -s ions.tpr -o ${name}_packed_ions.gro -p topol.top -pname NA -nname CL -neutral

echo '1 q'|gmx make_ndx -f ${name}_packed_ions.gro -o system.ndx

# Minimization 1
gmx grompp -f ${mdp_dir}/em.mdp -c ${name}_packed_ions.gro -p topol.top -n system.ndx -o em.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -v -deffnm em

# 100 ps NVT with protein chains frozen in place using freeze groups
gmx grompp -f ${mdp_dir}/nvt.mdp -c em.gro -r em.gro -p topol.top -n system.ndx -o nvt.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -deffnm nvt -nb gpu -pme gpu -bonded gpu

# Initial water molecules placed in the box using packmol now seep into the gaps between the protein grid
# Re-solvate the system to ensure proper solvation
gmx solvate -cp nvt -o nvt_solvated.gro -p topol.top

echo '1 q'|gmx make_ndx -f nvt_solvated.gro -o system.ndx

# Minimization 2 to relax any bad contacts after re-solvation
gmx grompp -f ${mdp_dir}/em.mdp -c nvt_solvated.gro -p topol.top -n system.ndx -o em2.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -v -deffnm em2

# 100 ps NVT after re-solvation with protein chains frozen in place using freeze groups
gmx grompp -f ${mdp_dir}/nvt.mdp -c em2.gro -r em2.gro -p topol.top -n system.ndx -o nvt_solvated.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -deffnm nvt_solvated -nb gpu -pme gpu -bonded gpu

# Minimization 3 just to be safe before proceeding to NPT
gmx grompp -f ${mdp_dir}/em.mdp -c nvt_solvated.gro -p topol.top -n system.ndx -o em3.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -v -deffnm em3

# 100 ns semi-isotropic NPT with 1 bar pressure coupling in the xy plane and 0 bar in the z direction
# Protein chains are no longer frozen in place
gmx grompp -f ${mdp_dir}/npt.mdp -c em3.gro -r em3.gro -p topol.top -n system.ndx -o npt_test.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -deffnm npt_test -nb gpu -pme gpu -bonded gpu

# 300 ns semi-isotropic NPT with 1 bar pressure coupling in the xy plane and 0 bar in the z direction
# Protein chains are no longer frozen in place
gmx grompp -f ${mdp_dir}/npt2.mdp -c npt_test.gro -r npt_test.gro -p topol.top -n system.ndx -o npt2.tpr -maxwarn 3
gmx mdrun -ntmpi 1 -deffnm npt2 -nb gpu -pme gpu -bonded gpu