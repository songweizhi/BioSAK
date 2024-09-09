library(ape)
library(GMD)
library(tools)
library(optparse)

# usgae
# Rscript ~/R_scripts/newick_tree/Tree2grouping.R -t human_gut_tree.newick

MOAR_LETTERS <- function(n=2) {
  n <- as.integer(n[1L])
  if(!is.finite(n) || n < 2)
    stop("'n' must be a length-1 integer >= 2")
  res <- vector("list", n)
  res[[1]] <- LETTERS
  for(i in 2:n)
    res[[i]] <- c(sapply(res[[i-1L]], function(y) paste0(y, LETTERS)))
  unlist(res)
}
group_id_letters = MOAR_LETTERS(3)

option_list = list(
  
  make_option(c("-t", "--tree"), 
              type="character", 
              default=NULL, 
              help="input file name", 
              metavar="character"));
  
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
wd = getwd()

tree_file_in = opt$tree
tree_file_in_name_no_extension = file_path_sans_ext(basename(tree_file_in))
grouping_file_out = paste(tree_file_in_name_no_extension, 'grouping.txt', sep = '_')
tree_txt_file_with_group = paste(tree_file_in_name_no_extension, 'with_group.txt', sep = '_')
tree_txt_file_only_group = paste(tree_file_in_name_no_extension, 'only_group.txt', sep = '_')
tree_plot_file = paste(tree_file_in_name_no_extension, '.jpg', sep = '_')
tree_plot_file_with_group = paste(tree_file_in_name_no_extension, 'with_group.jpg', sep = '_')
tree_plot_file_only_group = paste(tree_file_in_name_no_extension, 'only_group.jpg', sep = '_')

pwd_grouping_file_out = paste(wd, grouping_file_out, sep = '/')
pwd_tree_txt_file_with_group = paste(wd, tree_txt_file_with_group, sep = '/')
pwd_tree_txt_file_only_group = paste(wd, tree_txt_file_only_group, sep = '/')

pwd_tree_plot_file = paste(wd, tree_plot_file, sep = '/')
pwd_tree_plot_file_with_group = paste(wd, tree_plot_file_with_group, sep = '/')
pwd_tree_plot_file_only_group = paste(wd, tree_plot_file_only_group, sep = '/')

# read in tree file
SCG_tree = read.tree(tree_file_in)

# get distance matrix from input tree
distance_matrix = cophenetic(SCG_tree)

# sort the matrix
distance_matrix = distance_matrix[order(row.names(distance_matrix)),order(row.names(distance_matrix))]

# get the optimal k value
dist.obj <- dist(distance_matrix)
hclust.obj <- hclust(dist.obj)
css.obj <- css.hclust(dist.obj,hclust.obj)
elbow.obj <- elbow.batch(css.obj, ev.thres = 0.9)
#elbow.obj <- elbow.batch(css.obj)
k <- elbow.obj$k; cutree.obj <- cutree(hclust.obj,k=k)

# perform kmeans analysis of the matrix and get predicted clusters
clusters = kmeans(distance_matrix, k, nstart = 1000)$cluster

# turn predicted clusters to dataframe and write out 
clusters_df = data.frame(clusters)

# trun row name (bin id) to name column
clusters_df_row_name_in_column <- data.frame(names = row.names(clusters_df), clusters_df, row.names = NULL)

# sort dataframe according to group index
clusters_df_sorted = clusters_df_row_name_in_column[order(clusters_df_row_name_in_column$clusters),]

# add group information 
clusters_df_sorted$group = with(clusters_df_sorted, group_id_letters[clusters])

# extract need information for output
clusters_df_for_write = data.frame(clusters_df_sorted$group, clusters_df_sorted$names)

# write out into grouping .txt 
write.table(clusters_df_for_write, file = pwd_grouping_file_out, row.names=FALSE, col.names = FALSE, sep = ',', quote = FALSE)

# plot input tree
SCG_tree = read.tree(tree_file_in)
jpeg(pwd_tree_plot_file, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree, 'u', font = 1, cex = 0.5, label.offset = 0, lab4ut = 'axial')
dev.off()


# plot SCG tree with group
SCG_tree_with_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_with_group$tip.label)) {
  label_name = SCG_tree_with_group$tip.label[i]
  label_name_row_num = which(clusters_df_for_write$clusters_df_sorted.names == label_name)
  group_id = clusters_df_for_write$clusters_df_sorted.group[label_name_row_num]
  SCG_tree_with_group$tip.label[i] = paste(group_id, SCG_tree_with_group$tip.label[i], sep = '_')
  i = i + 1}
# plot tree with group
jpeg(pwd_tree_plot_file_with_group, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_with_group, 'u', font = 1, cex = 0.5, label.offset = 0, lab4ut = 'axial')
dev.off()
write.tree(SCG_tree_with_group, file=pwd_tree_txt_file_with_group)

# plot SCG tree only group
SCG_tree_only_group = read.tree(tree_file_in)
i = 1
for (i in 1:length(SCG_tree_only_group$tip.label)) {
  label_name = SCG_tree_only_group$tip.label[i]
  label_name_row_num = which(clusters_df_for_write$clusters_df_sorted.names == label_name)
  group_id = clusters_df_for_write$clusters_df_sorted.group[label_name_row_num]
  SCG_tree_only_group$tip.label[i] = paste(group_id)
  i = i + 1}
# plot tree only group
jpeg(pwd_tree_plot_file_only_group, width = 1600, height = 1600, units = "px", quality = 100, pointsize = 30)
plot.phylo(SCG_tree_only_group, 'u', font = 1, cex = 0.7, label.offset = 0, lab4ut = 'axial')
#plot.phylo(SCG_tree_only_group, font = 1, cex = 0.7, label.offset = 0.1, lab4ut = 'axial')

dev.off()
write.tree(SCG_tree_only_group, file=pwd_tree_txt_file_only_group)



