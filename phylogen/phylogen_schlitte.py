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
from dendropy import Tree, TreeList, DnaCharacterMatrix, DataSet, TaxonNamespace

#____________________________________________________________________________
# Read phylogenetic data
#____________________________________________________________________________

t = 'consensus.tree'
a = 'C1.P1.initial.fasta'

def itree(tree, alignment):
	itree.tree = dendropy.Tree.get(path=tree, schema="Newick")
	bipart = itree.tree.encode_bipartitions()
	itree.aa1 = dendropy.ProteinCharacterMatrix.get(file=open(alignment), schema="fasta")

for node in itree.tree:
	node.edge.bipartition.leafset_taxa(itree.taxa)
