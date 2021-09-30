
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)
}


# install packages if not
packages<-c("optparse", "ggplot2")
invisible(suppressMessages(check.packages(packages)))


############################################### Usage example ##############################################

# !!! COG annotation results for the HGTs must be in the last one row
# Rscript ~/Dropbox/R_scripts/for_plot/COG_summary_boxplot_last1row.R -i matrix.csv -o plot.png

############################################################################################################

option_list = list(
  
  make_option(c("-i", "--input"), 
              type="character", 
              help="input file name", 
              metavar="character"),
  
  make_option(c("-o", "--out"), 
              type="character", 
              help="output plot name", 
              metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


############################################### read in data ###############################################

my_data = read.csv(opt$input, row.names = 1)

# remove columns with all zero values
my_data = my_data[, colSums(my_data != 0) > 0]

# slice dataframe
row_num = dim(my_data)[1]
input_genomes_COG = my_data[1:(row_num-1),]
HGT_COG = my_data[row_num,]


######################################### get shape and color list #########################################

shape_list  = c()
color_list  = c()
for(each_col in names(my_data)){
  each_col_quatile = quantile(my_data[each_col][1:(row_num-1),])
  each_col_quatile_low = each_col_quatile[[2]]
  each_col_quatile_high = each_col_quatile[[4]]
  each_col_HGT = my_data[each_col][row_num,]
  
  current_shape = ''
  current_color = ''
  if ((each_col_HGT - each_col_quatile_high) > 0){
    current_shape = 24
    current_color = 'firebrick2'
  } else if ((each_col_HGT - each_col_quatile_low) < 0){
    current_shape = 25
    current_color = 'deepskyblue'
  } else {
    current_shape = 22
    current_color = 'black'
  }
  # add current shape and color to their list
  shape_list = c(shape_list, current_shape)
  color_list = c(color_list, current_color)
}

################################################# get plot #################################################

png(filename=opt$out, units="in", width=9, height=6, pointsize=10, res=300)

ggplot() + 
  geom_boxplot(data = stack(my_data[1:(row_num-1),]), aes(x = ind, y = values)) + 
  
  geom_point(data = data.frame(x = factor(colnames(my_data)), 
                               y = as.numeric(my_data[row_num,])),
             aes(x=x, y=y),
             color = color_list, fill = color_list, 
             shape = shape_list,
             stroke = 0.8, size = 2.5) + 
  
  theme_bw() + # white background
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black",size=6,angle=270,hjust=0,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=6,angle=0,hjust=0,vjust=0,face="plain"),
        axis.title.x = element_text(colour="black",size=6,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=6,angle=90,hjust=.5,vjust=0,face="plain"))+
  labs(x = "Category", y = 'Proportion') # label of x and y axis

invisible(dev.off())

