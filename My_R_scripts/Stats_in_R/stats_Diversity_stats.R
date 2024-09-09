
#################### Diversity Indices for Community Ecologists ####################

# load library
library(vegan)


?diversity


# define a dataframe 
species_abundance = c(2, 8, 1, 1, 3)
data_frame = data.frame(species_abundance)
row.names(data_frame) = c('apple', 'banana', 'orange', 'peach', 'pear')
data_frame


# Shannon diversity
diversity(data_frame, index = 'shannon')


# Simpson diversity
diversity(species_abundance, index = 'simpson')
diversity(data_frame, index = 'simpson')

library(diverse)
data_frame_t = t(data_frame)
data_frame_t
diversity(data_frame_t, type = "simpson")

