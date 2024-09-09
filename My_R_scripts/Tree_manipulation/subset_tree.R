# install phyloseq
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

library(phylobase)
library(phyloseq)
library(ape)



data(geospiza)
nodeLabels(geospiza) <- paste("N", nodeId(geospiza, "internal"), sep="")
geotree <- extractTree(geospiza)

class(geotree)
plot(geotree)

GTDB_tree = readNewick('/Users/songweizhi/Desktop/arc122_r86.1.tree')
SCG_tree = readNewick('/Users/songweizhi/Desktop/GoodBins_0.5_0.05_species_tree.newick')

class(SCG_tree)

plot(SCG_tree)
?readNewick

# tips <- c("difficilis", "fortis", "fuliginosa", "fusca", "olivacea", "pallida", "parvulus", "scandens")
tips <- c('psittacula', 'parvulus', 'pallida', 'pauper', 'fusca')
tree_subset = subset(geotree, tips.include=tips)
plot(tree_subset)


GTDB_tree = read.tree('/Users/songweizhi/Desktop/arc122_r86.1.tree')
class(GTDB_tree)

taxon_to_keep = c('GB_GCA_002725245.1', 'GB_GCA_002496725.1', 'GB_GCA_002687305.1')
GTDB_tree_subset = subset(GTDB_tree, tips.include=taxon_to_keep)


?subset
plot(GTDB_tree_subset)


GTDB_tree_subset2 = prune_samples(taxon_to_keep, GTDB_tree)

??prune_samples
