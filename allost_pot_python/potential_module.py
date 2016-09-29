#############################################################################
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
##############################################################################
##                  IMPORTS
##############################################################################
from __future__ import unicode_literals
from future.utils import bytes_to_native_str
from prody import *
from pylab import *
import argparse
import numpy as np
import numpy.random as npr
import subprocess
import itertools as it
from numpy import linalg as LA
import math as math
import array as ar
import matplotlib.pyplot as plt

##############################################################################
##               Allosteric potential calculation                           ##
##############################################################################

protein_filename ='protein'
ref_conf_filename ='ref_coord'
modes_free_filename = 'modesfree.dat'
modes_bound_filename = 'modesrestr.dat'

protein = np.loadtxt(protein_filename)
ref_conf = np.loadtxt(ref_conf_filename)
modes_free = np.loadtxt(modes_free_filename)
modes_bound = np.loadtxt(modes_bound_filename)
            
def AlloPotFree(protein,ref_conf,modes_free,modes_bound):    

    kT=0.6 
    mmode=16 
    
    AlloPot_free={i: np.zeros(mmode) for i in protein}
    AlloPot_bound={i: np.zeros(mmode) for i in protein}
    
    for a1, a2 in it.combinations(protein,2):
        a1 = int(a1)
        a2 = int(a2)
        dX = ref_conf[int(a1-1)]-ref_conf[int(a2-1)] 
        d0 = LA.norm(dX)
        if d0<11.0: # Filter for first neighbours
            for m in range(mmode):
              AlloPot_free[a1][m] += LA.norm(modes_free[m][a1]-modes_free[m][a2])**2 
              AlloPot_free[a2][m] += LA.norm(modes_free[m][a1]-modes_free[m][a2])**2
              AlloPot_bound[a1][m] += LA.norm(modes_bound[m][a1]-modes_bound[m][a2])**2
              AlloPot_bound[a2][m] += LA.norm(modes_bound[m][a1]-modes_bound[m][a2])**2

    deltaGatom={resi: 0.5*kT*sum(return_array(AlloPot_bound[resi],AlloPot_free[resi])) for resi in protein}
    return deltaGatom


def return_array(bound, free):
    result = []
    for i in range(len(bound)):
        result.append(math.log(bound[i]/float(free[i])))
        return result
    
    

deltaG = AlloPotFree(protein,ref_conf,modes_free,modes_bound) 

#write result to file
def writedic(dic):
    with open('deltaG.dat', 'w') as f:
        for key, value in dic.items():
            f.write('%s  %s\n' % (key, value))

writedic(deltaG)

#########################################################
## Plot result
#########################################################
def plotdG(n):
    plt.bar(range(len(n)), n.values(), align='center')
    plt.xticks(np.arange(min(n), max(n)+1, 45.0))
    plt.ylabel(r'$\Delta G (kcal/mol)$')
    plt.xlabel(r'$C_\alpha$ atoms')
    plt.grid()
    axes = plt.gca()
    axes.set_xlim([min(n)-10,max(n)+10])
    axes.set_ylim([-1.2,1.2])
    plt.savefig('deltaG_plot.pdf')

plotdG(deltaG)
