library(ape)
library(tools)
library(optparse)



# read in tree file
SCG_tree = read.tree('/Users/songweizhi/Desktop/species_tree_with_grouping_3.newick')

# get distance matrix from input tree
distance_matrix = cophenetic(SCG_tree)

distance_matrix

write.csv(distance_matrix, '/Users/songweizhi/Desktop/species_tree_with_grouping_3.matrix')
