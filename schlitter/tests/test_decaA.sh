#! /bin/sh
#valgrind ../src/schlitter \
../src/schlitter \
--pdb decaA_1_calpha.pdb \
--traj decaA_1.g96 \
--coarse \
|| exit 1

