library("vegan")
library("pheatmap")

rm(list = ls())
mydir = getwd()

########################################################################################

# file in
file_nm     = 'tree_subset_33_Nitrosopumilaceae_Matrix_distance.txt'
plot_height = 30
plot_width  = 30

# file out
pdf_out_BC  = 'tree_subset_33_Nitrosopumilaceae_Matrix_distance_Heatmap_BC_BC.pdf'
txt_out_BC  = 'tree_subset_33_Nitrosopumilaceae_Matrix_distance_Heatmap_BC_BC.txt'
pdf_out_Eu  = 'tree_subset_33_Nitrosopumilaceae_Matrix_distance_Heatmap_Eu_Eu.pdf'
txt_out_Eu  = 'tree_subset_33_Nitrosopumilaceae_Matrix_distance_Heatmap_Eu_Eu.txt'

########################################################################################

df_infile  = read.delim(paste(mydir, '/',file_nm, sep = ''), row.names=1, check.names=FALSE) # check.names=FALSE will not add X in front of a numeric ID.
heatmap_in = df_infile

##########################################################################################

# 用pheatmap函数画热图----Eulidean & Eulidean

# cube root transformation (euclidean)
drows_mod = dist(heatmap_in, method = "euclidean")
dcols_mod = dist(t(heatmap_in), method = "euclidean")

# generate new plot:
pdf(pdf_out_Eu, width = plot_width, height = plot_height)
pheatmap(heatmap_in, cluster_row = T, clustering_distance_rows = drows_mod, clustering_distance_cols = dcols_mod, clustering_method = "complete",
         cellwidth = 60, cellheight = 60, fontsize = 30, treeheight_row = 400, treeheight_col = 400)
dev.off()

pheatmap.dt <- pheatmap(heatmap_in, cluster_row = T, clustering_distance_rows = drows_mod, clustering_distance_cols = dcols_mod, clustering_method = "complete",
                        legend_breaks = 0:1, legend_labels = c("0","1"), fontsize = 3, treeheight_row = 100, treeheight_col = 100)

pheatmap.map.order <- (heatmap_in)[pheatmap.dt$tree_row$order,pheatmap.dt$tree_col$order]
write.table(x = pheatmap.map.order, file = txt_out_Eu, append = F, quote = F, sep = "\t", row.names = T, col.names = T)

########################################################################################

# 用pheatmap函数画热图----BC & BC

# cube root transformation (BC distance)
drows_mod_BC <- vegdist(heatmap_in, method="bray", binary=FALSE, diag=FALSE, upper=FALSE,na.rm = F)
dcols_mod_BC <- vegdist(t(heatmap_in), method="bray", binary=FALSE, diag=FALSE, upper=FALSE,na.rm = F)

# generate new plot:
pdf(pdf_out_BC, width = plot_width, height = plot_height)
pheatmap(heatmap_in, cluster_row = T, clustering_distance_rows = drows_mod_BC, clustering_distance_cols = dcols_mod_BC, clustering_method = "complete",
         cellwidth = 60, cellheight = 60, fontsize = 30, treeheight_row = 400, treeheight_col = 400)
dev.off()

pheatmap.dt.BC_BC <- pheatmap(heatmap_in, cluster_row = T, clustering_distance_rows = drows_mod_BC, clustering_distance_cols = dcols_mod_BC, clustering_method = "complete",
                              legend_breaks = 0:1, legend_labels = c("0","1"),fontsize = 5, treeheight_row = 100, treeheight_col = 100)

pheatmap.map.order.BC_BC <- heatmap_in[pheatmap.dt.BC_BC$tree_row$order,pheatmap.dt.BC_BC$tree_col$order]
write.table(x = pheatmap.map.order.BC_BC, file = txt_out_BC, append = F, quote = F, sep = "\t", row.names = T, col.names = T)
