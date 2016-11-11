#! /bin/sh
#valgrind ../src/schlitter \
../src/schlitter \
--pdb test_short_traj.pdb \
--traj test_short_traj.g96 \
--coarse \
|| exit 1

