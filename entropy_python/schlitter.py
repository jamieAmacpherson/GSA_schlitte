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
from MDAnalysis.analysis.align import *
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

parser.add_argument("-t", dest="dcdfile", required=True,
                    help="Free trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))


parser.add_argument("-s", dest="pdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))


args = parser.parse_args()

#____________________________________________________________________________
# remove rotational-translational motions from trajectory
#____________________________________________________________________________
def rmrt(topology, trajectory):
    ref = mda.Universe(topology)
    traj = mda.Universe(topology, trajectory)
    rms_fit_trj(traj, ref, filename='rmsfit_traj.dcd')

rmrt(args.pdbfile, args.dcdfile)

#____________________________________________________________________________
# covariance matrix
#____________________________________________________________________________
def covar(topology, trajectory):
    struct = parsePDB(topology)
    traj = Trajectory(trajectory)
    traj.link(struct)
    traj.setCoords(struct)
    traj.setAtoms(struct.calpha)
    ensemble = EDA('trajectory')
    ensemble.buildCovariance( traj )
    mat = ensemble.getCovariance()
    covar.mat = mat 
    np.savetxt('covarmat.dat', covar.mat )

covar(args.pdbfile, 'rmsfit_traj.dcd')
  
#____________________________________________________________________________
# mass matrix
#____________________________________________________________________________
def mass(topology, trajectory):
    system = mda.Universe(topology, trajectory)
    asel = system.select_atoms(' ( name CA ) ')
    am = asel.masses * 1.6605e-27 # convert from atomic mass to kg
    masses = np.repeat(am, 3)
    mass.matrix = np.identity(len(masses)) * masses

mass(args.pdbfile, args.dcdfile)

#____________________________________________________________________________
# entropy (reduced units)
#____________________________________________________________________________
def entropy(sigma, m):
    hbar = 6.6260693e-34 / (2 * np.pi)  # J*s
    k = 1.380658e-23  #J/K
    n = 6.0221367e23 # mol
    T = 300  # Kelvin
    be = (k * T * math.exp(2) / (hbar**2))
    print "units:", units
    print "logdet:", LA.slogdet((np.identity(len(m)) + (be * np.dot(sigma, m))))[1]
    
    entropy.Su = (k/2 * n * ((LA.slogdet((np.identity(len(m)) + (units * np.dot(sigma, m)))))[1]))
    print "S': ", entropy.Su, "J/(mol K)"
    
entropy(covar.mat, mass.matrix)






