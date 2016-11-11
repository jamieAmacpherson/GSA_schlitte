#!/bin/bash
if [ $# -ne 3 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.pdb)> time step (ps)"
	exit 1
fi

schlitter='/home/jamie/jm.software/GSA_schlitte/entropy_python/src/schlitter.py'



let k=$3
let nst=500000

while [ $k -le $nst ]
do



trjconv -f $1 -s $2 -b 0 -e $k -o trajout$k.xtc <<EOF
3
EOF

trjconv -f $1 -s $2 -b 0 -e 0 -o topol$k.pdb <<EOF
3
EOF

mkdir $k
cd $k
python ../../mdconvert.py ../trajout$k.xtc -o $k.dcd

 
python $schlitter -t $k.dcd -s ../topol$k.pdb 

rm ../trajout$k.xtc $k.dcd ../topol$k.pdb

cd ..
let k=k+$3
done
exit
