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
import fnmatch
import glob
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

# Link the phylogenetic tree to fasta sequences 
def itree(tree, alignment):
	itree.taxa = dendropy.TaxonNamespace()
	itree.tree = dendropy.Tree.get(path=tree, schema="Newick", taxon_namespace=itree.taxa)
	bipart = itree.tree.encode_bipartitions()

itree(t,a)

# Extract bipartitions
def getbipart(tree):
	getbipart.biparts=[]
	# extract bipartitions in the tree and append to empty list
	for node in itree.tree:
		getbipart.biparts.append(node.edge.bipartition.leafset_taxa(itree.taxa))
	# save the list of tree bipartitions to a txt file
	np.savetxt('bipartition.dat', getbipart.biparts, delimiter='', fmt="%s")
	os.system('chmod +wx format_newick.sh')
	os.system('./format_newick.sh')

getbipart(t)


def matchnewicks(alignment):
	# for each of the bipartitions import the sequence header,
	# match it to the aligned sequence and save that
	# sequence as a dataframe
	# Import the alignment file and save as a tmp object
	fastatmp = open(alignment, 'r').read().splitlines()
        # Remove double spaces from tmp alignment file
	fastatmp2 = [x.strip('  ') for x in fastatmp]
        # Replace underscore with single space in sequence headers
	fastatmp3 = [w.replace('_', ' ') for w in fastatmp2]
	# Remove '>' from the begining of each sequence header
	# in the alignment file.
	matchnewicks.fasta = [w.replace('>', '') for w in fastatmp3]	
	#
	matchnewicks.bpheaders = {}
	# Import bipartition sequence headers into a dictionary
	# maintaining the hierarchy within the tree topology (ie.
	# 1.seq corresponds to common ancestral sequence and n.seq
	# corresponds to the header sequence 
	for filename in glob.glob('*.inseq'):
		f_contents = open(filename, 'r').read().strip().splitlines()
		matchnewicks.bpheaders[filename] = f_contents
	#
	#  
		for seq, headers in matchnewicks.bpheaders.iteritems():
			matchnewicks.result=[]
			for header in headers:
				for i, element in enumerate(matchnewicks.fasta):
					if header in element:
						matchnewicks.result.append(matchnewicks.fasta[i+1])
				for filename in matchnewicks.bpheaders:
					matchnewicks.bpheaders[filename] = headers, matchnewicks.result

matchnewicks(a)

