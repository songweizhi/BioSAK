library(gplots)
library(ape)
#dat <- read.tree(file="your/newick/file")
dat <- read.tree(text="((A:4.2,B:4.2):3.1,C:7.3);")

#####################################################################################################
## All Pairwise Intersects: to identify and cluster overlaps among two to thousands of sample sets ##
#####################################################################################################
setlist <- lapply(1:6, function(x) sample(letters, x, replace=TRUE)); names(setlist) <- paste("S", seq(along=setlist), sep="")

setlist

# To work with the following function, the sample sets (here 6) need to be stored in a list object (here 'setlist')
# where each component is named by a unique identifier (here S1-S6).
setlist <- sapply(setlist, unique)

setlist

olMA <- sapply(names(setlist), function(x) sapply(names(setlist), function(y) sum(setlist[[x]] %in% setlist[[y]]))); olMA
# Creates an intersect matrix for all pairwise sample comparisons stored in 'setlist'
olMA

typeof(olMA)
?heatmap.2
#heatmap.2(olMA, trace="none", col=colorpanel(40, "darkred", "orange", "yellow"), Colv='Rowv', Rowv="none", dendrogram="none")
heatmap.2(olMA, trace="none", col=colorpanel(40, "darkred", "orange", "yellow"), Colv='Rowv', Rowv="none")
heatmap.2(olMA, trace="none", col=colorpanel(40, "darkred", "orange", "yellow"))

# Hierarchical clustering of the rows and columns of the intersect matrix 'olMA'. The result is plotted as heatmap
