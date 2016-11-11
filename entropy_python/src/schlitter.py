#____________________________________________________________________________
# Entropy calculation for an MD simulation using the Schlitter equation
#
# Copyright 2016 King's College London and the Authors
# 
# Author: Jamie A. Macpherson, Jens Kleinjung and Franca Fraternali
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
import matplotlib.pyplot as plt
import csv
import mpmath as mp

#____________________________________________________________________________
# Parse commandline arguments
#____________________________________________________________________________
# check if input file exists
# there are two inputs: trajectory (.dcd) and topology (.pdb)
# if either of those inputs are not supplied, or if the user doesn't invoke the
# help flag the program will display an error message.
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

# command line argument parser
parser = argparse.ArgumentParser(description='Calculate entropy of MD trajectory')

# the first argument is the trajectory file (.dcd) supplied after the -t flag
# the trajectory file is saved as an object with the variable args.dcdfile
parser.add_argument("-t", dest="dcdfile", required=True,
                    help="Free trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))

# the second argument is the topology file (.pdb) supplied after the -s flag
# this is saved an an obect with the variable args.pdbfile
parser.add_argument("-s", dest="pdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))

# the arguments are parsed 
args = parser.parse_args()

#____________________________________________________________________________
# remove rotational-translational motions from trajectory
#____________________________________________________________________________
# to calculate the configurational entropy, we first remove the rotational-translational
# motions from the trajectory. This is done by fitting the trajectory to the reference
# structure (ie. the pdbfile). 

def rmrt(topology, trajectory):
    # define the reference structure as the topology file
    ref = mda.Universe(topology)
    
    # define the trajectory as the .dcd file and link it to the topology symbolically
    traj = mda.Universe(topology, trajectory)
    
    # fit the trajectory to the topology, removing rotational-translational motions
    # and save the resulting trajectory.
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
    
    # the covariance matrix is calculated in square angstrom (1 angstrom^2 = 0.01 nm^2)
    covar.mat = mat * 0.01	# convert non-mass weighted covariance matrix from a^2 to nm^2  
    np.savetxt('covarmat_nonMW.dat', covar.mat)
    
    # Mass-weight the covariance matrix, this involves a simple factoring of each of the elements
    # by the atomic mass unit of carbon and hydrogen. The resulting units are covariance in: AMU*nm^2
    covar.matMW = (covar.mat * 13.01864)	# weight the atoms by atomic mass of H + C (united atom force field)
    np.savetxt('covarmat_MW.dat', covar.matMW)
    
    covar.matMWSI = (covar.matMW * 1.66054e-45)	# convert covariance to SI units from U*nm^2 to kg*m^2
    np.savetxt('covarmat_MW_SI.dat', covar.matMWSI)

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
# entropy
#____________________________________________________________________________
# calculate entropy using the Schlitter equation

def entropy(sigma): 
    
    # define constants all in SI units
    hbar = 6.6260693e-34 / (2 * np.pi)  # J*s
    k = 1.380658e-23  #J/K
    n = 6.0221367e23 # mol
    T = 300  # Kelvin
    
    # the Schlitter equation is S' = 0.5k_{B} \sum_{i=1}^{N=3} (1 + \frac{k_{B} T e^{2}}{\hbar^{2}} * \langle q \rangle) 
    # determine the product of the prefactors first 
    be = (k * T * math.exp(2) / (hbar**2))
  
    # diagonalize the mass-weighted matrix and calculate its eigenvalues  
    eigenvals, eigenvects = LA.eig(sigma)
    
    # for each eigenvector, calculate the inner product of the Schlitter equation and sum over all eigenvectors
    for key in eigenvals:
	    deter = []
	    dd = mp.log(1 + (be * key))	# natural log (using mpmath module - enables long float operations) 
	    deter.append(dd)
    logdeter = np.sum(deter)
    
    # multiply the log inner product by the Boltzmann constant  
    entropy.S = 0.5 * k * n * logdeter
    
    # write the calculated entropy to an output file
    f = open( 'entropy.dat', 'w' )
    f.write( str(entropy.S) )
    f.close()

    print "S': ", entropy.S, "J/(mol K)"
    
    # Measure entropy summing over different number of eigenvalues
    entropy.moderange={ }
    for rmode in range(len(eigenvals)):
	    tmp_arr = []
	    rmdd = mp.log(1 + be * sum(eigenvals[0:rmode]))
	    rmS = 0.5 * k * n * rmdd
	    tmp_arr.append(rmS)
            entropy.moderange[rmode] = tmp_arr
    
    with open('entropy_moderange.csv', 'w') as f:
	    c = csv.writer(f)
	    for key, value in entropy.moderange.items():
		    c.writerow([key] + value)
   
    # plot eigenvectors of the covariance matrix
    cumeigval = np.cumsum(eigenvals)
    plt.plot(cumeigval, 'b .', cumeigval, 'r--')
    plt.ylabel('Sum of eigenvalues')
    plt.xlabel('Sum of eigenvectors')
    plt.grid()
    plt.savefig('eigenval_spectrum.pdf')
    
    
entropy(covar.matMWSI)
