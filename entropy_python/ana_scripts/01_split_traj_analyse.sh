#!/bin/bash
if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: (time step (ps)) (last time step (ps))"
	exit 1
fi

plot='/home/jamie/jm.software/GSA_schlitte/entropy_python/ana_scripts/split_traj_plot.R'

let k=$1
let nst=$2

while [ $k -le $nst ]
do

cp $k/entropy.dat entropy$k".dat"

awk 'FNR==1{print ""}1' entropy$k.dat >> time_entropies.dat

let k=k+$1
done

rm entropy**.dat

R CMD BATCH $plot

exit
