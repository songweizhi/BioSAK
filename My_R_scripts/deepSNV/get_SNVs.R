
# example command
# Rscript ~/R_scripts/deepSNV/get_SNVs.R -q D1.bam -c D0.bam -i test_ctg -s 1 -e 300
# Rscript ~/R_scripts/deepSNV/get_SNVs.R -q 2D9.bam -c D2D0.bam -i D2_c -s 1 -e 3000000

############################## install library ##############################

# # install GenomeInfoDbData first
# source("https://bioconductor.org/biocLite.R")
# biocLite("GenomeInfoDbData")
# 
# # update packages
# install.packages('matrixStats')
# install.packages('VGAM')
# update.packages('matrixStats')
# update.packages('VGAM')
# 
# # install deepSNV (http://master.bioconductor.org/packages/devel/bioc/html/deepSNV.html)
# # try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("deepSNV")

# load library
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(deepSNV)))

################################## optparse #################################

option_list = list(
  
  make_option(c("-q", "--query_bam"), 
              type="character", 
              default=NULL, 
              help="query bam file", 
              metavar="character"),
  
  make_option(c("-c", "--control_bam"), 
              type="character", 
              default=NULL, 
              help="control bam file", 
              metavar="character"),
  
  make_option(c("-i", "--ctg_id"), 
              type="character", 
              default=NULL, 
              help="contig id", 
              metavar="character"),
  
  make_option(c("-s", "--start_pos"), 
              type="integer", 
              default=NULL, 
              help="start position", 
              metavar="number"),
  
  make_option(c("-e", "--end_pos"), 
              type="integer", 
              default=NULL, 
              help="end position", 
              metavar="number"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
query_bam = opt$query_bam
control_bam = opt$control_bam
ctg_id = opt$ctg_id
start_pos = as.integer(opt$start_pos)
end_pos = as.integer(opt$end_pos)

# define file names
query_bam_path = dirname(query_bam)
query_bam_filename_no_extension = file_path_sans_ext(basename(query_bam))
control_bam_filename_no_extension = file_path_sans_ext(basename(control_bam))
SNV_output_jpg_file = paste(query_bam_filename_no_extension, 'vs', control_bam_filename_no_extension, ctg_id, start_pos, end_pos, 'SNVs.jpg', sep = '_')
SNV_output_txt_file = paste(query_bam_filename_no_extension, 'vs', control_bam_filename_no_extension, ctg_id, start_pos, end_pos, 'SNVs.txt', sep = '_')

################################## get SNV ##################################

# set regions
regions <- data.frame(chr=ctg_id, start = start_pos, stop = end_pos)
SNV_simulated <- deepSNV(test = query_bam, control = control_bam, regions=regions)
# plot(SNV_simulated)
# png(filename=SNV_output_jpg_file, units="in", width=30, height=30, pointsize=12, res=300)
png(filename=SNV_output_jpg_file)
plot(SNV_simulated)
invisible(dev.off())

# get summary
# SNVs <- summary(SNV_simulated, sig.level=0.05, adjust.method="BH")
SNVs <- summary(SNV_simulated)

# write out
write.table(SNVs, file = SNV_output_txt_file, row.names=FALSE, col.names = TRUE, sep = ',', quote = FALSE)

