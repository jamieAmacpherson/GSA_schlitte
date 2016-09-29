#____________________________________________________________________________
# Entropy calculation for an MD simulation using the Schlitter equation
#
# Copyright 2016 King's College London and the Authors
# 
# Author: Jamie A. Macpherson
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
from prody import *
from pylab import *
import MDAnalysis as mda
import argparse
import os.path
import sys
import numpy as np
import subprocess
import itertools as it
from numpy import linalg as LA
import math as math
import array as ar

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

# covariance matrix
def covar(topology, trajectory):
    struct = parsePDB(topology)
    traj = Trajectory(trajectory)
    traj.link(struct)
    traj.setCoords(struct)
    traj.setAtoms(struct.calpha)
    ensemble = EDA('trajectory')
    ensemble.buildCovariance( traj )
    covar.covarmat = ensemble.getCovariance()
    
covar(args.pdbfile, args.dcdfile)

# mass matrix
def mass(topology, trajectory):
    system = mda.Universe(topology, trajectory)
    asel = system.select_atoms(' ( name CA ) ')
    am = asel.masses
    masses = np.repeat(am, 3)
    mass.matrix = np.identity(len(masses)) * masses

mass(args.pdbfile, args.dcdfile)

# entropy (reduced units)
def entropy_ru(sigma, m):
    entropy_ru.S = LA.logsdet(np.dot(sigma, m))[1]

entropy_ru(covar.covarmat, mass.matrix)

def entropy(sigma, m):
    h = 6.6260693e-34 / (2 * np.pi)  # J*s
    k = 0.0000083144621  # J/(mol⋅K)
    T = 300  # Kelvin
    units = (k * T * math.exp(2) / (h**2))
    entropy.Su = (k/2 * ((LA.logsdet((1 + units * np.dot(sigma, m))))[1]))
    print entropy.Su, "J/(mol K)"
    
    with open("config_entropy.txt","w") as f:
        for line in entropy.Su[0]:
            strs=" ".join(str(x) for x in line)
            f.write(strs+"\n")

entropy(covar.covarmat, mass.matrix)






