set sel [atomselect top all]
set M [transvecinv {27.3680 105.1850 -51.8760}] 
$sel move $M
set M [transaxis y -90]
$sel move $M
$sel writepdb ../output/elp_car9_aligned.pdb
quit