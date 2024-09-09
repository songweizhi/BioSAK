
# usage example
# cd /Users/songweizhi/Dropbox/Research/Flow_cell_datasets/deepSNV_output_summary
# Rscript /Users/songweizhi/PycharmProjects/FlowCellBiofilm/R_scripts/perform_MDS_and_PERMANOVA.R -s 210 -m deepSNV_output_summary_210_frequency.txt -f deepSNV_output_summary_210_factor.txt -v euclidean -p deepSNV_output_summary_210_frequency_MDS.png


# check.packages function: install and load multiple R packages.
# Check to see if packages are installed. Install them if they are not, then load them into the R session.
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = 1)
}


# library vegan for metaMDS
# library ecodist for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds
packages = c("optparse", "vegan", "ecodist")
check.packages(packages)


options(warn=-1)

option_list = list(
  
  make_option(c("-s", "--strain"), 
              type="character", 
              help="strain, 210/D2", 
              metavar="character"),
  
  make_option(c("-m", "--matrix"), 
              type="character", 
              help="input matrix", 
              metavar="character"),
  
  make_option(c("-f", "--factor"), 
              type="character", 
              help="input matrix factor", 
              metavar="character"),
  
  make_option(c("-v", "--vegdist"), 
              type="character", 
              help="vegdist method, euclidean/jaccard/bray", 
              metavar="character"),
  
  make_option(c("-p", "--plot"), 
              type="character", 
              help="plot file name", 
              metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


####################################### Main code #######################################

# Import data
snv_summary = read.delim(opt$matrix, row.names=1)
snv_factor = read.delim(opt$factor, row.names=1)

# check samples and variables, # 638 variables (snv), 36 samples
#dim(snv_summary)

# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)
# calculate relative abundance
# snv_summary_t_norm = snv_summary_t / rowSums(snv_summary_t) * 100
## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = opt$vegdist)
# plot hclust
# snv_summary_t_ja_hclus = hclust(snv_summary_t_ja, method = "average") # Cluster the distances
# par(pty = 's') # square plot
# plot(snv_summary_t_ja_hclus)
# plot nMDS
# non-metric multi-dimensional scaling (nMDS), which is an ordination method.
# A low stress value (<0.2) and a linear pattern in a stressplot suggest a good mapping of your data.
#snv_summary_t_ja_mds = metaMDS(snv_summary_t_ja, autotransform = F, trace = F, trymax=50)
snv_summary_t_ja_mds = metaMDS(comm = snv_summary_t_ja)
#snv_summary_t_ja_mds = nmds(snv_summary_t_ja)

# plot stress
# stressplot(snv_summary_t_ja_mds)


####################################### Get MDS plot #######################################

# get plot_shape_list, Mono210: square; MonoD2: triangle; Coculture: circle
plot_shape_list  = c()
for (eachA in snv_factor$Species){
  if (eachA == 'Mono210'){
    current_shape = 15
  } else if (eachA == 'MonoD2'){
    current_shape = 17
  } else if (eachA == 'Coculture'){
    current_shape = 19
  }
  plot_shape_list = c(plot_shape_list, current_shape)
}

# plot
png(filename=opt$plot, units="in", width=10, height=10, pointsize=12, res=300)

par(mar=c(5,4,2,7))
plot(snv_summary_t_ja_mds$points, col=snv_factor$Label, pch=plot_shape_list)
text(snv_summary_t_ja_mds$points, labels=snv_factor$Time, cex= 0.7, pos=3)

# get legend_shape_list
legend_shape_list = c()
legend_shape_list_210 = c(19, 19, 19, 15, 15, 15)
legend_shape_list_D2 = c(19, 19, 19, 17, 17, 17)
if (opt$strain == '210'){
  legend_shape_list = legend_shape_list_210
}else if (opt$strain == 'D2'){
  legend_shape_list = legend_shape_list_D2
}

# add legend
legend("topleft", inset=c(1,0), xpd=TRUE, bty="n", legend = levels(snv_factor$Label),pch=legend_shape_list, col=as.numeric(as.factor(levels(snv_factor$Label))))


############################## Perform PERMANOVA analysis ##############################

# PERMANOVA analysis (http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/adonis.html)
adonis2(snv_summary_t_ja ~ snv_factor$Species, snv_summary, method = vegdist_method)

# clear variables/objects
rm(list=ls())

