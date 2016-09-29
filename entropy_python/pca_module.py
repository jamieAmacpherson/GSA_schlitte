#____________________________________________________________________________
# Allosteric potential calculation from residue first neighbors d0 = 1.1 nm
#
# Copyright 2016 King's College London and the Authors
# 
# Authors: James Macpherson
# Contributors: Enrico Guarnera, Anna Ladach
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
from prody import *
from pylab import *
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

# redirect stdout to file
sys.stdout = open('pca.out', 'w')

# check if input file exists
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

parser = argparse.ArgumentParser(description='Allosteric potential calculation from residue first neighbors d0 = 1.1 nm')

parser.add_argument("-ft", dest="freedcdfile", required=True,
                    help="Free trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))


parser.add_argument("-fs", dest="freepdbfile", required=True,
                    help="Free structure file (format: .pdb)",
                    type=lambda x: is_valid_file(x))


parser.add_argument("-rt", dest="restdcdfile", required=True,
                    help="Restrained trajectory file (format: .dcd)",
                    type=lambda x: is_valid_file(x))


args = parser.parse_args()

##

#TO DO: parse filenames as strings instead of objects. 
#
#____________________________________________________________________________

def trajpca_superimposed(freepdbfile, freedcdfile, restdcdfile):
    freestruct = parsePDB(freepdbfile)
    freedcd = Trajectory(freedcdfile)
    nframes = freedcd.numCoordsets()
    freedcd.addFile(restdcdfile)
    freedcd.link(freestruct)
    freedcd.setCoords(freestruct)
    freedcd.setAtoms(freestruct.calpha)

    freetraj = freedcd[:nframes]
    freetraj.superpose()

    restraj = freedcd[nframes:]
    restraj.superpose()

    pcafree = PCA('PCA analysis free')
    freecovmat1 = pcafree.buildCovariance(freetraj)
    pcafree.calcModes()
    writeNMD('freepca.nmd', pcafree[:20], freestruct.select('calpha'))

    pcarest = PCA('PCA analysis restrained')
    restrcovmat1 = pcarest.buildCovariance(restraj)
    pcarest.calcModes()
    writeNMD('restpca.nmd', pcarest[:20], freestruct.select('calpha'))

trajpca_superimposed(args.freepdbfile, args.freedcdfile, args.restdcdfile)

#TO DO: superimpose into the same coordinate system (ie. into same PCA space)
#
###################################################################################################
##Read in the pca files and manipulate them as inputs for allosteric potential calculation
#using bash commands to reorder the output of the previous section for modes, these can then be
#directly read into the allosteric potential calculation

def extractmodes(nmafree, nmarestr):
    #free modes
    extractmodes.free = []
    for line in nmafree:
        if line[:4] == 'mode':
            extractmodes.free.append(line[7:])
    extractmodes.freenp = np.asarray(extractmodes.free)
    np.savetxt('modesfree.dat', extractmodes.freenp, delimiter='', fmt="%s")
    #restr modes
    extractmodes.restr = []
    for line in nmarestr:
        if line[:4] == 'mode':
            extractmodes.restr.append(line[7:])
    extractmodes.restrnp = np.asarray(extractmodes.restr)
    np.savetxt('modesrestr.dat', extractmodes.restrnp, delimiter='', fmt="%s")


freepcafile = open("freepca.nmd", "r")    
restrpcafile = open("restpca.nmd", "r")
extractmodes(freepcafile, restrpcafile)


#extract protein residues
def extractcoord(pdb):
    extractcoord.resid = []
    extractcoord.xyz = []
    for line in pdb:
        if line[:4] == 'ATOM' and line[12:16] == " CA ":
            extractcoord.resid.append(line[23:26])
            extractcoord.xyz.append(line[32:54])
    extractcoord.xyz=np.asarray(zip(*[iter(extractcoord.xyz)]*1))
    extractcoord.resid=np.asarray(zip(*[iter(extractcoord.resid)]*1))
    np.savetxt('ref_coord', extractcoord.xyz, delimiter='', fmt="%s")
    np.savetxt('protein', extractcoord.resid, delimiter='', fmt="%s")

pdb = open(args.freepdbfile, "r")

extractcoord(pdb)



