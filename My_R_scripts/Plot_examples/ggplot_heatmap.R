library(ggplot2)
library(ggcorrplot)
library(reshape)

mydata <- mtcars[, c(1,3,4,5,6,7)]

head(mydata)

cormat <- round(cor(mydata),2)

head(cormat)
typeof(cormat)

setwd('/Users/songweizhi/Desktop/heatmap_test')

mydata = read.csv('ANI.txt', header = TRUE, sep=",",row.names = 1)
mydata

p.mat <- cor_pmat(mydata)
ggcorrplot(mydata)



#mydata_matrix <- as.matrix(mydata)
#class(mydata_matrix)
#heatmap(mydata_matrix)

# Generate decending values from 100 to 1 to simulate retention over time
rcohorts <- replicate(15, sort(runif(15, 1, 100), T))

# Make a triangle
rcohorts[lower.tri(rcohorts)] <- NA

# Convert to a data frame, and add tenure labels
rcohorts <- as.data.frame(rcohorts)
rcohorts$tenure <- seq(0,14)

# Reshape to suit ggplot, remove NAs, and sort the labels
rcohorts <- na.omit(melt(rcohorts, 'tenure', variable_name='cohort'))
rcohorts$cohort <- factor(rcohorts$cohort, levels=rev(levels(rcohorts$cohort)))
#The data shows a particular value for a given cohort across a number of weeks.
head(rcohorts)


ggplot(rcohorts, aes(cohort, tenure)) +
  ggtitle('Retention by cohort') +
  theme_bw() +
  xlab('Cohort') +
  ylab('Tenure (weeks)') +
  geom_tile(aes(fill = value), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab') +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'))

