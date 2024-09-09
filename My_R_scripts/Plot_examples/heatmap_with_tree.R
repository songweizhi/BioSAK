library(gplots)   # contains the heatmap.2 package
library(car)   

States
dim(States)

scaled <- scale(States[,-1]) # scale all but the first column to make information comparable
heatmap.2(scaled, # specify the (scaled) data to be used in the heatmap
          cexRow=0.5, cexCol=0.95, # decrease font size of row/column labels
          scale="none", # we have already scaled the data
          trace="none") # cleaner heatmap



