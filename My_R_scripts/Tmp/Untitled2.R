library(ape)
library(tools)
library(optparse)

tree_file_in = '/Users/songweizhi/Desktop/bac120_r83.tree.txt'
tree_plot_out = '/Users/songweizhi/Desktop/bac120_r83.tree.jpg'

  
# plot SCG tree only group
SCG_tree = read.tree(tree_file_in)

# plot tree only group
jpeg(tree_plot_out, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree, 'u', font = 1, cex = 0.7, label.offset = 0.1, lab4ut = 'axial')
dev.off()

