library(vegan) # for metaMDS
library(ecodist) # for nmds,  https://www.rdocumentation.org/packages/ecodist/versions/2.0.1/topics/nmds

# set working directory
setwd("/Users/songweizhi/Desktop/000")

vegdist_method = 'euclidean' # euclidean or jaccard or bray, manhattan
# Dissimilarity index, partial match to "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis".
# Jaccard is useful for presence/absence analyses (composition only)

# get input file name
strain = '210' # 210 or D2
snv_summary_file = paste('SNV_QC_ncd_even_flk_depth_', strain, '_matrix_no_plasmid2.txt', sep = '')
snv_factor_file = paste('Stats_', strain, '_factor.txt', sep = '')


##############################################

work_dir          = '/Users/songweizhi/Desktop/Biofilm2021/OneStep_MinBoth_10_MinEach_1_StrandBias_10_DepthDiff_30'
snv_summary_file  = 'SNV_QC_ncd_even_flk_depth_210_matrix.txt'
snv_factor_file   = 'Stats_210_factor.txt'
strain            = '210' # 210 or D2
vegdist_method    = 'euclidean'


##############################################

setwd(work_dir)

# plot out
plot_out = paste(snv_summary_file, '.', vegdist_method, '.pdf', sep = '')

# Import data
snv_summary = read.delim(snv_summary_file, row.names=1)
snv_factor = read.delim(snv_factor_file, row.names=1)

# transpose dataframe as vegan wants variables (SNVs) as columns
snv_summary_t = t(snv_summary)

## Calculate dissimilarity matrix, type ?vegdist for other options (e.g. Jaccard is useful for presence/absence analyses (composition only)).
snv_summary_t_ja = vegdist(snv_summary_t, method = vegdist_method)
snv_summary_t_ja_mds = metaMDS(comm = snv_summary_t_ja)

####################################### Get MDS plot #######################################

plot_label_list = c()
plot_color_list = c()
plot_shape_list = c()
plot_label_list_uniq = c()
plot_color_list_uniq = c()
plot_shape_list_uniq = c()
for (each_point in row.names(snv_summary_t_ja_mds$points)){
  sample_code = substring(each_point, 1, 2)
  sample_code = strsplit(each_point, 'D')[[1]][1]

  if (sample_code == 'X4'){
    current_label = 'Coculture_A'
    current_color = 'darkorange'
    current_shape = 19
  } else if (sample_code == 'X8'){
    current_label = 'Coculture_B'
    current_color = 'darkorchid1'
    current_shape = 19
  } else if (sample_code == 'X12'){
    current_label = 'Coculture_C'
    current_color = 'red'
    current_shape = 19
  } else if (sample_code == 'X1'){
    current_label = 'Mono210_A'
    current_color = 'blue1'
    current_shape = 15
  } else if (sample_code == 'X5'){
    current_label = 'Mono210_B'
    current_color = 'deepskyblue'
    current_shape = 15
  } else if (sample_code == 'X9'){
    current_label = 'Mono210_C'
    current_color = 'chartreuse1'
    current_shape = 15
  } else if (sample_code == 'X2'){
    current_label = 'MonoD2_A'
    current_color = 'blue1'
    current_shape = 17
  } else if (sample_code == 'X6'){
    current_label = 'MonoD2_B'
    current_color = 'deepskyblue'
    current_shape = 17
  } else if (sample_code == 'X10'){
    current_label = 'MonoD2_C'
    current_color = 'chartreuse1'
    current_shape = 17
  }
  
  plot_label_list = c(plot_label_list, current_label)
  plot_color_list = c(plot_color_list, current_color)
  plot_shape_list = c(plot_shape_list, current_shape)
  
  if (current_label %in% plot_label_list_uniq == FALSE){
    plot_label_list_uniq = c(plot_label_list_uniq, current_label)
    plot_color_list_uniq = c(plot_color_list_uniq, current_color)
    plot_shape_list_uniq = c(plot_shape_list_uniq, current_shape)
  }
}


# plot
pdf(plot_out, width=5, height=5, pointsize=5)
par(mar=c(5,4,2,7))
plot(snv_summary_t_ja_mds$points, col=plot_color_list, pch=plot_shape_list)
text(snv_summary_t_ja_mds$points, labels=snv_factor$Time, cex= 0.7, pos=3)
legend("topleft", inset=c(1,0), xpd=TRUE, bty="n", legend = plot_label_list_uniq, pch=plot_shape_list_uniq, col=plot_color_list_uniq)
invisible(dev.off())
