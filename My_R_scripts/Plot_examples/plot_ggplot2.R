library(ggplot2) #load package
library(datasets) #load R datasets
library(gcookbook) #load gcookbook datasets 

heightweight #dataset

# Plots:
ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 1) + geom_point(size = 2) #shape: geom_point(shape = 2), fill:geom_point(size = 1)
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, color = sex)) + geom_point() #p77

# set 
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, color = sex)) + geom_point() + 
  scale_shape_manual(values = c(1,3))
  scale_color_brewer(palette = 'set1') #p77
  
#use slightly larger points and use a shape scale with custom values
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point(size =3) + scale_shape_manual(values = c(1,2))
