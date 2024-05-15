from ete3 import Tree
import sys

# Load the PHYLIP tree
phylip_tree = Tree(sys.argv[1], format=1)

# Write the tree to Newick format
phylip_tree.write(format=0, outfile=sys.argv[2])

