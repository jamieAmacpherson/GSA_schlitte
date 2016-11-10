#!/bin/bash
if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.pdb)>"
	exit 1
fi

schlitter='/home/jamie/jm.software/GSA_schlitte/entropy_python/src/schlitter.py'



let k=8000
let nst=300000

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

let k=k+8000
done
exit

for d in */ ; do
	cat $d/entropy.dat >> time_entropy.dat
done

