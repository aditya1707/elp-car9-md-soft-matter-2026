#!/bin/bash

name="elp_car9_aligned"
#create a grid of ELPs

file="packmol_input.inp"
touch ${file}

x=10
y=10
counter=1
echo "tolerance 2.0" >> ${file}
echo " " >> ${file}
echo "filetype pdb" >> ${file}
echo " " >> ${file}
#echo "output elp_r4qr12q_aligned_packed.pdb" >> ${file}
echo "output ${name}_packed.pdb" >> ${file}
echo " " >> ${file}
echo "pbc 50.0 50.0 120.0" >> ${file}
for i in {1..3}
do
    for j in {1..3}
    do
	echo $x $y $counter
	echo " " >> ${file}
#	echo "structure elp_r4qr12q_aligned_final.pdb" >> ${file}
	echo "structure ${name}_final.pdb" >> ${file}
	echo "  number 1" >> ${file}
	#echo "  atoms 435" >> ${file}
	#echo "      center" >> ${file}
	echo "      fixed ${x} ${y} 0. 0. 0. 0." >> ${file}
	#echo "      inside box 0. 0. 0. 50. 50. 110." >> ${file}
	echo "  end atoms" >> ${file}
	echo "end structure">>${file}
	echo " " >> ${file}
	y=$(( y + 15 ))
	counter=$(( counter + 1 ))
    done
    y=10
    x=$(( x + 15 ))
done

echo "structure water.pdb" >> ${file}
echo "number 2500" >> ${file}
echo "inside box 0. 0. 73. 50. 50. 120." >> ${file}
echo "end atoms" >> ${file}
echo "end structure" >> ${file}

