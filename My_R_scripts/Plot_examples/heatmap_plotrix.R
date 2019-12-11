
library(plotrix)

input_file = "/Users/songweizhi/Desktop/heatmap_plotrix_example_data.csv"

chromosome = read.csv(input_file, header = TRUE)
chromosome_matrix = as.matrix(chromosome[,-1])
rownames(chromosome_matrix) = chromosome[,1]

row_number = dim(chromosome_matrix)[1]
column_number = dim(chromosome_matrix)[2]

color2D.matplot(chromosome_matrix, 
                xlab="Genome", 
                ylab="COG category", 
                show.legend=T, 
                border=NA, 
                color.spec="rgb", 
                axes=FALSE,
                cs1=c(0,1), cs2=c(0,0), cs3=c(0,0) # control the color, see below a few examples
                )

axis(1,at=0.5:(column_number-0.5),labels=colnames(chromosome_matrix)) 
axis(2,at=(row_number - 0.5):0.5, labels=rownames(chromosome_matrix)) 

rownames(chromosome_matrix)
# black and red: cs1=c(0,1), cs2=c(0,0), cs3=c(0,0)
# archaea yellow: cs1=c(0,1), cs2=c(0,1), cs3=c(0,0)
# gamma purple: cs1=c(0,1), cs2=c(0,0), cs3=c(0,1)
# gamma red+yellow: cs1=c(1,1), cs2=c(1,0), cs3=c(0,0)
# c(0,1),c(1,0.6,0.4,0.3,0),c(0.1,0.6)

# clear variables/objects
rm(list=ls())
