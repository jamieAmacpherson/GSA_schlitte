#!/bin/bash
if [ $# -ne 5 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.tpr)> <time step (ps)> <final time step (ps)> <index file (.ndx)>"
	exit 1
fi

set -e

extr="/home/macphej/jm.software/development/entropy/GSA_schlitte/entropy_python/ana_scripts/entropy_extract.py"

let k=$3
let nst=$4

while [ $k -le $nst ]
do

gmx trjconv -f $1 -s $2 -b 0 -e $k -o trajout$k.xtc -n $5 -fit rot+trans <<EOF
3
3
EOF

gmx trjconv -f $1 -s $2 -b 0 -e 0 -o topol$k.pdb -n $5 <<EOF
3
3
EOF

mkdir $k
cd $k

gmx covar -f ../trajout$k.xtc -s ../topol$k.pdb <<EOF
3
3
EOF

gmx anaeig -entropy > entropycal.dat
 
python $extr

rm ../trajout$k.xtc ../topol$k.pdb 

cd ..
let k=k+$3
done
exit
