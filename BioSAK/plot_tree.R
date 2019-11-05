
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Usage example
packages<-c("ape", "tools", "optparse")
invisible(suppressMessages(check.packages(packages)))


# usgae
# Rscript ~/R_scripts/newick_tree/Tree2grouping.R -t human_gut_tree.newick

option_list = list(
  
  make_option(c("-t", "--tree"), 
              type="character", 
              default=NULL, 
              help="input file name", 
              metavar="character"));


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


tree_file_in = opt$tree
tree_file_in_name_no_extension = file_path_sans_ext(basename(tree_file_in))
pwd_tree_plot_file = paste(tree_file_in_name_no_extension, 'pdf', sep = '.')

# read in tree file
SCG_tree = read.tree(tree_file_in)

# plot input tree
SCG_tree = read.tree(tree_file_in)

SCG_tree_leaf_number = SCG_tree$tip.label
tree_height = length(SCG_tree_leaf_number)*15

pdf(pwd_tree_plot_file, width = 1600, height = tree_height)
#png(pwd_tree_plot_file, width = 1000, height = tree_height, pointsize = 30)
#jpeg(pwd_tree_plot_file, width = 1600, height = tree_height, units = "px", quality = 100, pointsize = 30)
#plot.phylo(SCG_tree, 'u', font = 1, cex = 1, label.offset = 0.01, lab4ut = 'axial')
plot.phylo(SCG_tree, type = 'phylogram', show.tip.label = TRUE, font = 1, cex = 0.5, no.margin = FALSE, align.tip.label = TRUE, plot = TRUE)  # Type: phylogram,cladogram, fan, unrooted, radial.
add.scale.bar(x= -1, y = 0)
dev.off()


