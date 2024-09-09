



# http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual?tmpl=%2Fsystem%2Fapp%2Ftemplates%2Fprint%2F&showPrintDialog=1
# 
# #############################################################################################
# ## All Possible Intersects: for comprehensive overlap analyses of 2-20 or more sample sets ##
# #############################################################################################
# source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") 
# # To compute all possible intersects among more than two samples at the same time - such as common in AB, ABC, 
# # ABCD and so on - the overLapper function can be used. Note: processing more than 20 sets with this function can take 
# # a long time, because the complexity of all possible intersects increases with the number of sample sets n according 
# # to this relationship: (2^n) - 1 (includes self comparisons). This means there are 63 possible intersects for n=6, 1,023 
# # for n=10, 32,767 for n=15 and 33,554,431 for n=25. With the current implementation, the computation time is about 0.1 
# # second for n=10 and 2.5 seconds for n=15.
# setlist <- lapply(1:6, function(x) sample(letters, 20)); names(setlist) <- LETTERS[seq(along=setlist)] 
# # To work with the overLapper function, the sample sets need to be stored in a list object where the different 
# # components are named by unique identifiers: A, B, C, ..., Z.
# OLlist <- overLapper(setlist=setlist, complexity=1:length(setlist), sep="-", type="intersects"); OLlist; names(OLlist) 
# # Computes all possible intersects for the samples stored in 'setlist'. The results are returned as a list where each 
# # overlap component is labeled by the corresponding sample names. For instance, the name 'A-B-C' indicates that the 
# # entries are common in samples A, B and C. The seperator used for naming the intersect samples can be specified under 
# # the 'sep' argument. The complexity level range to consider can be controlled with the 'complexity' argument, which 
# # can have values from 2 to the total number of samples.
# OLlist[[2]]; OLlist[[3]] # Returns the corresponding intersect matrix and complexity levels.
# OLexport <- as.matrix(unlist(sapply(OLlist[[4]], paste, collapse=" ")))
# write.table(OLexport, file="test.xls", col.names=F, quote=F, sep="\t") # Exports intersect data in tabular format to a file.
# counts <- sapply(OLlist[[4]], length) # Counts the number of elements in each intersect component.
# tapply(counts, OLlist[[3]], function(x) rev(sort(x))) 
# # Sorts the overlap results within each complexity level by their size. This allows to identify the sample set 
# # combinations with the largest intersect within each complexity level.
# x11(height=12, width=8); olBarplot(OLlist=OLlist, horiz=T, las=1, cex.names=0.6, main="Intersect Bar Plot") 
# # Plots the counts as bar diagram. More details on this function are provided in the Venn diagram section. 
 




#####################################################################################################
## All Pairwise Intersects: to identify and cluster overlaps among two to thousands of sample sets ##
#####################################################################################################
setlist <- lapply(11:30, function(x) sample(letters, x, replace=TRUE)); names(setlist) <- paste("S", seq(along=setlist), sep="")
# To work with the following function, the sample sets (here 20) need to be stored in a list object (here 'setlist')
# where each component is named by a unique identifier (here S1-S20).
setlist <- sapply(setlist, unique)
olMA <- sapply(names(setlist), function(x) sapply(names(setlist), function(y) sum(setlist[[x]] %in% setlist[[y]]))); olMA
# Creates an intersect matrix for all pairwise sample comparisons stored in 'setlist'. This approach scales well up to
# 3several thousands of sample sets.
library("gplots")
#heatmap.2(olMA, trace="none", Colv="none", Rowv="none", dendrogram="none", col=colorpanel(40, "darkred", "orange", "yellow")) # Plots intersect matrix as heatmap.
heatmap.2(olMA, trace="none", col=colorpanel(40, "darkred", "orange", "yellow"))
# Hierarchical clustering of the rows and columns of the intersect matrix 'olMA'. The result is plotted as heatmap
# with two identical dendrograms representing the outcome of the hierarchical clustering. The latter is internally
# performed by calls of heatmap.2() to the functions dist() and hclust() using their default settings: euclidean
# distances and complete linkage. Note: the distance matrix used for clustering in this examples is based on the
# row-to-row (column-to-column) similarities in the olMA object. The next example shows how one can use the olMA
# object directly as distance matrix for clustering after transforming the intersect counts into similarity measures,
# such as Jaccard or Rand indices.
# diffMA1 <- sapply(names(setlist), function(x) sapply(names(setlist), function(y) sum(!setlist[[x]] %in% setlist[[y]])))
# diffMA2 <- sapply(names(setlist), function(x) sapply(names(setlist), function(y) sum(!setlist[[y]] %in% setlist[[x]])))
# jaccardMA <- olMA/(diffMA1+diffMA2+olMA)
# hr <- hclust(as.dist(1-jaccardMA), method = "complete", members=NULL)
# heatmap.2(jaccardMA, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hr), trace="none", col=colorpanel(40, "darkred", "orange", "yellow"))
# # Uses the olMA for clustering after transfroming its intersect counts into Jaccard (Tanimoto) similarity indices.
# # This transformation can give reasonable results for sample sets with large size differences. Many alternative
# # similarity measures for set comparisons could be considered here.
# setlist[["S1"]][setlist[["S1"]] %in% setlist[["S2"]]]
# # Example for returning the intersect entries for one pairwise comparison, here S1 and S2.
# name1 <- matrix(colnames(olMA), length(olMA[,1]), length(olMA[,1]), byrow=T)
# name2 <- matrix(colnames(olMA), length(olMA[,1]), length(olMA[,1]), byrow=F)
# mynames <- paste(name1, name2, sep="_"); mynames <- lapply(strsplit(mynames, "_"), sort)
# mynames <- sapply(mynames, paste, collapse="_"); olV <- as.vector(olMA)
# names(olV) <- mynames; olV <- olV[!duplicated(names(olV))]; sort(olV)
# # Converts the intersect matrix 'olMA' into a named vector sorted by intersect size. This allows to identify the
# # combinations of pairwise sample comparisons that have the largest or smallest intersects.
# 



#' ######################################################################################################
#' ## Presence-absence matrices: to indentify memberships of items across large numbers of sample sets ##
#' ######################################################################################################
#' setlist <- lapply(11:30, function(x) sample(letters, x, replace=TRUE)); names(setlist) <- paste("S", seq(along=setlist), sep="") 
#' # Generates a sample data set.
#' source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R") 
#' # Imports the required overLapper() function.
#' paMA <- overLapper(setlist=setlist, type="intersects", complexity=2)[[2]]
#' paMA[names(rev(sort(rowSums(paMA)))),] 
#' # Creates a presence-absence matrix of all items across their sample sets. In this matrix the presence information is 
#' # indicated by ones and its rows are sorted by the presence frequencies. Alternative versions of the present-absent 
#' @ matrix can be returned by setting the 'type' argument to 1, 2 or 3.
#' library("gplots"); heatmap.2(paMA, trace="none", Colv="none", Rowv="none", dendrogram="none", col=c("white", "gray")) 
#' # Plots the present-absent matrix as heatmap.
#' rowSums(paMA) # Returns the presence frequencies of the items in the sample sets.
#' sapply(rownames(paMA), function(x) colnames(paMA)[paMA[x, ]==1]) 
#' # Returns for each item the names of sample sets where it is present. The opposite result for absent calls can be 
#' # returned by changing '==1' to '==0'.








