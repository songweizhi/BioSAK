library(vegan)
library(ggplot2)

annotation_wd = '/Users/songweizhi/Desktop/Kelp_coassembly/annotation'
sample_metadata = paste(annotation_wd, 'sample_metadata.txt', sep = '/')

annotation_file = paste(annotation_wd, 'df_dbcan_norm.txt', sep = '/')
#annotation_file = paste(annotation_wd, 'df_COG_id_norm.txt', sep = '/')   # id, cate
#annotation_file = paste(annotation_wd, 'df_KEGG_B_norm.txt', sep = '/')   # B, C, D



dbcan_annotation <- read.csv(file = annotation_file, header = TRUE, sep = '\t')

sample_metadata = read.csv(file = sample_metadata, header = TRUE, sep = '\t')
dbcan_annotation_values = dbcan_annotation[2:ncol(dbcan_annotation)]

dbcan.mds <- metaMDS(comm = dbcan_annotation_values, distance = "bray", trace = FALSE, autotransform = TRUE)

MDS_xy <- data.frame(dbcan.mds$points)
MDS_xy$Location <- sample_metadata$Location
MDS_xy$Replicate <- sample_metadata$Replicate
MDS_xy$Month <- sample_metadata$Month
ggplot(MDS_xy, aes(MDS1, MDS2, color = Location, shape = Replicate)) + geom_point()
