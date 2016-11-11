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

mkdir $k

trjconv -f $1 -s $2 -b 0 -e $k -o $k/trajout$k.xtc <<EOF
3
EOF

trjconv -f $1 -s $2 -b 0 -e 0 -o $k/topol$k.pdb <<EOF
3
EOF


python ../mdconvert.py $k/trajout$k.xtc -o $k/$k.dcd

cd $k 
python $schlitter -t $k.dcd -s topol$k.pdb 
cd ..
rm $k/trajout$k.xtc

let k=k+$3
done
exit

