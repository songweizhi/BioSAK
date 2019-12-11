
#################### Diversity Indices for Community Ecologists ####################

# load library
library(vegan)

?procrustes


# define a dataframe 
library(ape)
library(vegan)

TREE1 = read.tree('/Users/songweizhi/Desktop/species_tree_1.0.newick')
TREE2 = read.tree('/Users/songweizhi/Desktop/species_tree_0.4.newick')

D1 = cophenetic(TREE1)
D1 = D1[order(row.names(D1)),order(row.names(D1))]
D2 = cophenetic(TREE2)
D2 = D2[order(row.names(D2)),order(row.names(D2))]

procrustes_summary = procrustes(D1, D2)
summary(procrustes_summary)

plot(procrustes_summary)
plot(procrustes_summary, kind=2)
residuals(procrustes_summary)




