#!/bin/bash
if [ $# -ne 2 ]
then
	echo "Incorrect number of arguments..."
	echo "Usage: split_resid.sh <trjectory (.xtc)> <topology (.pdb)>"
	exit 1
fi
python = '/usr/bin/python'
schlitter = '/Users/Jamie/jm.software/development/entropy/GSA_schlitte/entropy_python/src/schlitter.py'

let k=1
let nst=2

while [ $k -le $nst ]
do

mkdir $k
gmx make_ndx -f $2 -o $k/$k".ndx" <<EOF
r $k
q
EOF

gmx trjconv -f $1 -s $2 -n $k/$k".ndx" -o $k/$k.xtc <<EOF
14
EOF

mdconvert $k/$k.xtc -o $k/$k.dcd
cd $k $python $schlitter $k.dcd $2 
cd ..
rm $k/$k.xtc

let k=k+1
done
exit
