#!/bin/bash
gmx editconf -f npt2.gro -o npt2.pdb
plumed driver --plumed ${i}.dat --mf_xtc npt2.xtc --pdb npt2.pdb