library(ape)
library(vegan)

TREE1 = read.tree('/Users/songweizhi/Desktop/species_tree_denovo.newick')
TREE2 = read.tree('/Users/songweizhi/Desktop/species_tree_subset.newick')

D1 = cophenetic(TREE1)
D1 = D1[order(row.names(D1)),order(row.names(D1))]
D2 = cophenetic(TREE2)
D2 = D2[order(row.names(D2)),order(row.names(D2))]

mantel(xdis = D1, ydis = D2, permutations = 999)


