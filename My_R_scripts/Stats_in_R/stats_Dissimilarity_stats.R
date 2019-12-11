
#################### Dissimilarity Indices for Community Ecologists ####################

# load library
library(vegan)


# define a dataframe 
gold_fish = c(6, 10) 
guppies = c(7, 0) 
rainbow_fish = c(4, 6) 
data_frame = data.frame(gold_fish, guppies, rainbow_fish)
row.names(data_frame) = c('Tank1', 'Tank2')
data_frame

# Bray-Curtis dissimilarity
vegdist(data_matrix, method="bray")

