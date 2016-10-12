#____________________________________________________________________________
# Manipulating phylogenetic tree
#
# copyright 2016 king's college london and the authors
# 
# author: Jamie A. Macpherson, james.macpherson@crick.ac.uk
# 
# this library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with GSAtools. If not, see <http://www.gnu.org/licenses/>.
#____________________________________________________________________________
# Imports
#____________________________________________________________________________
import dendropy
import argparse
import os.path
import sys
import os
import numpy as np
#____________________________________________________________________________
# Parse commandline arguments
#____________________________________________________________________________

# check if input file exists
def is_valid_file(arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist! Use the --help flag for input options." % arg)
    else:
	return arg

parser = argparse.ArgumentParser(description='Model phylogenetic tree bipartitions')

parser.add_argument("-t", dest="treefile", required=True,
                    help="Phylogenetic tree file (Newick format)",
                    type=lambda x: is_valid_file(x))


parser.add_argument("-a", dest="alignmentfile", required=True,
                    help="Alignment file (Fasta format)",
                    type=lambda x: is_valid_file(x))


args = parser.parse_args()

#____________________________________________________________________________
# Read phylogenetic data
#____________________________________________________________________________

t = 'consensus.tree'
a = 'C1.P1.initial.fasta'

def itree(tree, alignment):
	itree.taxa = dendropy.TaxonNamespace()
	itree.tree = dendropy.Tree.get(path=tree, schema="Newick", taxon_namespace=itree.taxa)
	bipart = itree.tree.encode_bipartitions()

itree(t,a)

def getbipart(tree):
	getbipart.biparts=[]
	for node in itree.tree:
		getbipart.biparts.append(node.edge.bipartition.leafset_taxa(itree.taxa))
	np.savetxt('bipartition.dat', getbipart.biparts, delimiter='', fmt="%s")
	os.system('./format_newick.sh')
	#os.system('sed -n -e '/KPYR HUMAN/{p;n;}' bipartion.dat > human_bipart.dat')
	#os.system('awk -F\' '{ for (i=3; i<=NF; i+=2) $i=""} 1' human_bipart.dat > human_newicks1.dat')
	#os.system('cut -d " " -f 4- human_newicks1.dat > human_newicks.dat')

getbipart(t)

def matchnewicks(newicks):
	with open(newicks) as f:
		matchnewicks.seqhead = f.read().splitlines()
	for line in matchnewicks.seqhead:
		matchnewicks.seqhead[line].split()

matchnewicks('human_newicks.dat')

def twobytwo(t):
	it = iter(t)
	for x in it:
		yield x, next(it)

newheaders = dict(twobytwo(matchnewicks.seqhead))



