#!/bin/bash
if [ $# -ne 3 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.pdb)> time step (ps)"
	exit 1
fi

let k=$3
let nst=500000

while [ $k -le $nst ]
do

cp $k/entropy.dat entropy$k".dat"
awk 'FNR==1{print ""}1' entropy**.dat > time_entropies.dat

let k=k+$3
done

Rscript split_traj_plot.R 

exit
