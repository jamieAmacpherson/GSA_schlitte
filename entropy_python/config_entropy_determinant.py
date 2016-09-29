#____________________________________________________________________________
# Allosteric potential calculation from residue first neighbors d0 = 1.1 nm
#
# Copyright 2016 King's College London and the Authors
# 
# Author: James Macpherson
#
# GSAtools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with GSAtools. If not, see <http://www.gnu.org/licenses/>.
#____________________________________________________________________________
# Imports
#____________________________________________________________________________
import MDAnalysis as mda
from numpy import linalg as LA
import numpy as np
import math
import argparse
import os.path
import sys

#____________________________________________________________________________
# Parse commandline arguments
#____________________________________________________________________________
# check if input file exists
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

parser = argparse.ArgumentParser(description='Calculate entropy of MD trajectory')

parser.add_argument("-ft", dest="dcdfile", required=True,
                    help="Free trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))


parser.add_argument("-fs", dest="pdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))


args = parser.parse_args()

#____________________________________________________________________________
#  Import trajectory
#____________________________________________________________________________
topology = args.pdbfile
trajectory = args.dcdfile
system = mda.Universe(topology, trajectory)
asel = system.select_atoms(' ( name CA ) ')

#____________________________________________________________________________
#  making the mass matrix
#____________________________________________________________________________

kg = asel.masses * (1.6605388e-27)
masses = np.repeat(kg, 3) 
mass_matrix = np.identity(len(masses)) * masses
print mass_matrix
#____________________________________________________________________________
#  Preparing to read the CA position at every 5 steps in the traj
#____________________________________________________________________________

skip = 500
num_ts = system.trajectory.n_frames / skip
num_coor = len(asel) * 3
ca_pos = system.trajectory.timeseries(asel, skip=skip, format='fac')

#____________________________________________________________________________
# converting angstroms to meters and merging the xyz of timeseries
#____________________________________________________________________________

ca = (1e-10) * (np.reshape(ca_pos, (num_ts, -1)))
#print "ca", shape(ca)

#____________________________________________________________________________
#  making the covariance matrix
#____________________________________________________________________________

ca_avg = np.average(ca)
#print "ca_av", shape(ca_avg)
ca2 = ca - ca_avg[np.newaxis]
#print "ca2", shape(ca2)
ca_cov = np.zeros((num_coor, num_coor), dtype=float)
for ts in ca2:
    ca_cov += np.outer(ts, ts)
ca_cov /= num_ts
print ca_cov, np.shape(ca_cov)
print "mass_matrix", np.shape(mass_matrix)

#____________________________________________________________________________
#  calculating the entropy
#____________________________________________________________________________
hplanck_bar = 6.6260755e-34 / (2 * np.pi)  # J*s
k = (1.380658e-23)  # J/K
Avogadro = 6.0221367e23  # /mol
T = 300  # Kelvin
term = np.float64((k * T * math.exp(2) / pow(hplanck_bar, 2)))
print "term =", term
determ = LA.det((term * np.dot(ca_cov, mass_matrix)) + np.identity(len(mass_matrix)))
print "det = ", determ
S_ho = k / 2 * Avogadro * math.log(determ)
S_ho_kcal = S_ho * 0.000239006
print "S'=", S_ho_kcal, " kcal/(mol K)"


