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
		if node.edge.bipartition.leafset_as_newick_string(itree.taxa).find("KPYR") != -1:
			print node.edge.bipartition.leafset_as_newick_string(itree.taxa)
		
		getbipart.biparts.append(node.edge.bipartition.leafset_as_newick_string(itree.taxa)) 
	np.savetxt('bipartitions.dat', getbipart.biparts, delimiter='', fmt="%s")
	getbipart.matching = [s for s in getbipart.biparts if 'KPYR' in s]
	
getbipart(t)


