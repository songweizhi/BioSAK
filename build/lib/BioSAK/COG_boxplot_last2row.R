
check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = 1)
}


# install packages if not
packages<-c("optparse", "ggplot2")
invisible(suppressMessages(check.packages(packages)))


##################################### Usage example #####################################

# Rscript ~/Rscripts/COG_summary_boxplot.R -i COG_summary_47bins_in_percent.csv -o HGT_COG.png
# COG annotation results for the HGTs must be in the last two rows

#########################################################################################


get_shape_list_blast = function(csv_file){
  # read in dataframe from csv file
  my_data = read.csv(csv_file, row.names = 1)
  
  # get the number of rows and columns
  row_num = dim(my_data)[1]
  col_num = dim(my_data)[2]
  
  # slice dataframe
  input_genomes_COG = my_data[1:(row_num-2),]
  HGT_COG = my_data[row_num-1,]
  
  # get color list
  shape_list  = c()
  for(each_col in names(my_data)){
    #each_col_mean = mean(my_data[each_col][1:(row_num-2),])
    each_col_quatile = quantile(my_data[each_col][1:(row_num-2),])
    each_col_quatile_low = each_col_quatile[[2]]
    each_col_quatile_high = each_col_quatile[[4]]
      
    each_col_HGT = my_data[each_col][row_num-1,]
    #compare = each_col_HGT - each_col_mean
    current_shape = ''
    if ((each_col_HGT - each_col_quatile_high) > 0){
      current_shape = 24
    } else if ((each_col_HGT - each_col_quatile_low) < 0){
      current_shape = 25
    } else {
      current_shape = 22
    }
    shape_list = c(shape_list, current_shape)
  }
  return(shape_list)
}


get_shape_list_tree = function(csv_file){
  # read in dataframe from csv file
  my_data = read.csv(csv_file, row.names = 1)
  
  # get the number of rows and columns
  row_num = dim(my_data)[1]
  col_num = dim(my_data)[2]
  
  # slice dataframe
  input_genomes_COG = my_data[1:(row_num-2),]
  HGT_COG = my_data[row_num,]
  
  # get color list
  shape_list  = c()
  for(each_col in names(my_data)){
    each_col_mean = mean(my_data[each_col][1:(row_num-2),])
    each_col_quatile = quantile(my_data[each_col][1:(row_num-2),])
    each_col_quatile_low = each_col_quatile[[2]]
    each_col_quatile_high = each_col_quatile[[4]]
    
    each_col_HGT = my_data[each_col][row_num,]
    #compare = each_col_HGT - each_col_mean
    current_shape = ''
    if ((each_col_HGT - each_col_quatile_high) > 0){
      current_shape = 24
    } else if ((each_col_HGT - each_col_quatile_low) < 0){
      current_shape = 25
    } else {
      current_shape = 22
    }
    shape_list = c(shape_list, current_shape)
  }
  return(shape_list)
}


get_color_list_blast = function(csv_file){
  # read in dataframe from csv file
  my_data = read.csv(csv_file, row.names = 1)
  
  # get the number of rows and columns
  row_num = dim(my_data)[1]
  col_num = dim(my_data)[2]
  
  # slice dataframe
  input_genomes_COG = my_data[1:(row_num-2),]
  HGT_COG = my_data[row_num-1,]
  
  # get color list
  shape_list  = c()
  for(each_col in names(my_data)){
    #each_col_mean = mean(my_data[each_col][1:(row_num-2),])
    each_col_quatile = quantile(my_data[each_col][1:(row_num-2),])
    each_col_quatile_low = each_col_quatile[[2]]
    each_col_quatile_high = each_col_quatile[[4]]
    
    each_col_HGT = my_data[each_col][row_num-1,]
    #compare = each_col_HGT - each_col_mean
    current_shape = ''
    if ((each_col_HGT - each_col_quatile_high) > 0){
      current_shape = 'black'
    } else if ((each_col_HGT - each_col_quatile_low) < 0){
      current_shape = 'black'
    } else {
      current_shape = 'black'
    }
    shape_list = c(shape_list, current_shape)
  }
  return(shape_list)
}


get_color_list_tree = function(csv_file){
  # read in dataframe from csv file
  my_data = read.csv(csv_file, row.names = 1)
  
  # get the number of rows and columns
  row_num = dim(my_data)[1]
  col_num = dim(my_data)[2]
  
  # slice dataframe
  input_genomes_COG = my_data[1:(row_num-2),]
  HGT_COG = my_data[row_num,]
  
  # get color list
  shape_list  = c()
  for(each_col in names(my_data)){
    each_col_mean = mean(my_data[each_col][1:(row_num-2),])
    each_col_quatile = quantile(my_data[each_col][1:(row_num-2),])
    each_col_quatile_low = each_col_quatile[[2]]
    each_col_quatile_high = each_col_quatile[[4]]
    
    each_col_HGT = my_data[each_col][row_num,]
    #compare = each_col_HGT - each_col_mean
    current_shape = ''
    if ((each_col_HGT - each_col_quatile_high) > 0){
      current_shape = 'black'
    } else if ((each_col_HGT - each_col_quatile_low) < 0){
      current_shape = 'black'
    } else {
      current_shape = 'black'
    }
    shape_list = c(shape_list, current_shape)
  }
  return(shape_list)
}


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

csv_file = opt$input
my_data = read.csv(csv_file, row.names = 1)
row_num = dim(my_data)[1]
shape_list_tree = get_shape_list_tree(csv_file)
shape_list_blast = get_shape_list_blast(csv_file)
color_list_tree = get_color_list_tree(csv_file)
color_list_blast = get_color_list_blast(csv_file)

# plot
ggplot() + 
  geom_boxplot(data = stack(my_data[1:(row_num-2),]), aes(x = ind, y = values)) + 
  
  geom_point(data = data.frame(x = factor(colnames(my_data)), 
                               y = as.numeric(my_data[row_num-1,])),
             aes(x=x, y=y),
             color = color_list_blast, 
             shape = shape_list_blast,
             fill = color_list_blast, 
             stroke = 0.8,
             size = 2.5) + 
  
  geom_point(data = data.frame(x = factor(colnames(my_data)), 
                               y = as.numeric(my_data[row_num,])),
             aes(x=x, y=y),
             color = color_list_tree, 
             shape = shape_list_tree,
             stroke = 0.8,
             size = 2.5) + 

  theme_bw() + # white background
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="black",size=12,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.y = element_text(colour="black",size=12,angle=90,hjust=.5,vjust=0,face="plain"))+
  labs(x = "COG category", y = 'Proportion') # label of x and y axis

ggsave(opt$out, width = 10, height = 5)

