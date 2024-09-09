# library(optparse)
# library(ggplot2)
# library(plotrix)
# library(plsgenomics)
# library(ape)

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = 1)
}

# Usage example
packages<-c("optparse", "ggplot2", "plotrix", "plsgenomics", "ape")
check.packages(packages)


options(warn=-1)

option_list = list(
  
  make_option(c("-t", "--tree"), 
              type="character", 
              help="input GTDB tree subset", 
              metavar="character"),
  
  make_option(c("-m", "--matrix"), 
              type="character", 
              help="HGT matrix", 
              metavar="character"),
  
  make_option(c("-p", "--plot"), 
              type="character", 
              help="output plot", 
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

SCG_tree = read.tree(opt$tree)
HGT_df = read.csv(opt$matrix, sep=",", row.names = 1, header = TRUE) 
plot_file = opt$plot

# group level
# SCG_tree = read.tree("/Users/songweizhi/Desktop/subset_tree_c.newick")
# HGT_df = read.csv("/Users/songweizhi/Desktop/plot_c/HGT_matrix_group_level_norm.csv", sep=",", row.names = 1, header = TRUE) 
# plot_file = '/Users/songweizhi/Desktop/plot_c/HGT_matrix_group_level_norm.png'


# convert df into matrix
HGT_mat = as.matrix(HGT_df)
HGT_mat
# transpose dataframe
#HGT_mat = t(HGT_mat)
SCG_tree$tip
# re-order matrix labels according to SCG tree
HGT_mat_adj <- HGT_mat[SCG_tree$tip, SCG_tree$tip]


png(filename=plot_file, units="in", width=15, height=15, pointsize=12, res=300)

# set the position of x and y labels
xl <- c(0.5, ncol(HGT_mat_adj)+0.5)
yl <- c(0.5, nrow(HGT_mat_adj)+0.5)

# set layout
layout(matrix(c(0,1,0,2,3,4,0,5,0),nrow=3, byrow=TRUE),width=c(1,3,1), height=c(1,3,1))
#layout(matrix(c(0,0,1,0,0,0,2,0,3,4,5,6,0,0,7,0),nrow=4, byrow=TRUE),width=c(1,0.5,3,1), height=c(1,0.5,3,1))

# plot tree on top
par(mar=rep(0,4))
plot(SCG_tree, direction="downwards", show.tip.label=FALSE, xaxs="i", x.lim=xl)

# # group ID label on top
# par(mar=rep(0,4))
# #plot(NA, axes=FALSE, ylim=c(2,2), xlim=xl,  xaxs="i") #
# plot(NA, axes=FALSE, ylim=c(0,2), xlim=xl,  xaxs="i") #
# text(1:ncol(HGT_mat_adj),rep(2,ncol(HGT_mat_adj)), SCG_tree$tip, pos=4, srt = -90)

# plot tree on left
par(mar=rep(0,4))
plot(SCG_tree, direction="rightwards", show.tip.label=FALSE, yaxs="i", y.lim=rev(yl))

# # group ID label
# par(mar=rep(0,4))
# plot(NA, axes=FALSE, xlim=c(0,2), ylim=yl, yaxs="i") # 
# text(rep(0,nrow(HGT_mat_adj)),1:nrow(HGT_mat_adj),rev(SCG_tree$tip), pos=4)
 
     
# plot matrix
par(mar=rep(0,4), xpd=TRUE)
#image((1:nrow(HGT_mat))-0.5, (1:ncol(HGT_mat))-0.5, HGT_mat, axes=FALSE) # xaxs="i", yaxs="i"
#matrix.heatmap(HGT_mat)
#color2D.matplot(HGT_mat_adj, show.legend=F, border=NA, show.values=3, color.spec="rgb", axes=FALSE, cs1=c(0,1), cs2=c(1,0), cs3=c(0,0))
#color2D.matplot(HGT_mat_adj, show.legend=F, border=NA, show.values=3, color.spec="rgb", axes=FALSE)
color2D.matplot(HGT_mat_adj, show.legend=F, border=NA, show.values=0, color.spec="rgb", axes=FALSE)

# plot Y-axis labels
par(mar=rep(0,4))
plot(NA, axes=FALSE, xlim=c(0,2), ylim=yl, yaxs="i") # 
text(rep(0,nrow(HGT_mat_adj)),1:nrow(HGT_mat_adj),rev(SCG_tree$tip), pos=4)

# plot X-axis label
par(mar=rep(0,4))
#plot(NA, axes=FALSE, ylim=c(2,2), xlim=xl,  xaxs="i") #
plot(NA, axes=FALSE, ylim=c(0,2), xlim=xl,  xaxs="i") #
text(1:ncol(HGT_mat_adj),rep(2,ncol(HGT_mat_adj)), SCG_tree$tip, pos=4, srt = -90)

invisible(dev.off())
#rm(list=ls())

