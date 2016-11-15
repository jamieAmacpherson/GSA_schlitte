#!/bin/bash
if [ $# -ne 5 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.pdb)> <time step (ps)> <final time step (ps)> <index file (.ndx)>"
	exit 1
fi

schlitter='/Users/Jamie/jm.software/GSA_schlitte/entropy_python/src/schlitter.py'
#mdconvert='/Users/Jamie/jm.software/GSA_schlitte/entropy_python/ana_scripts/mdconvert.py'


let k=$3
let nst=$4

while [ $k -le $nst ]
do



gmx trjconv -f $1 -s $2 -b 0 -e $k -o trajout$k.xtc -n $5 <<EOF
3
EOF

gmx trjconv -f $1 -s $2 -b 0 -e 0 -o topol$k.pdb -n $5 <<EOF
3
EOF

mkdir $k
cd $k
python mdconvert ../trajout$k.xtc -o $k.dcd

 
python $schlitter -t $k.dcd -s ../topol$k.pdb 

rm ../trajout$k.xtc *.dcd ../topol$k.pdb

cd ..
let k=k+$3
done
exit
