# Paired t-tests
# http://environmentalcomputing.net/paired-t-tests/

# Paired t-tests are used to compare the means of two groups of measurements when individual objects are measured twice.

# For example, to contrast the photosynthetic performance of ten plants in two environments in a greenhouse (shady and sunny), 
# we could measure performance in each individual plant twice, once in the shade and once in the sun. 

# The 20 individual measurements are not independent of each other because we would expect the pair of measurements taken from 
# the same individual to be more similar to each other than if randomly sampled from all available plants. We are thus unable 
# to use an independent samples t-test – that test would be appropriate if each plant was used only once, with some plants
# measured in the shady treatment and different set of plants measured in the sunny treatment.

# Pairing data like this is usually done to reduce the likely variation among measurements with the aim of better detecting 
# differences between groups. In this example, the difference between the two measures of photosynthetic performance on a given 
# plant should reflect mostly the effect of sunlight, while in an independent samples design, the difference between a plant in 
# the shade and another plant in the sun will reflect both differences in the effects of sunlight and differences between the 
# individual plants.


############################################################
#                 Running the analysis                     #
############################################################

Greenhouse_data = '/Users/songweizhi/PycharmProjects/BioSAK/My_R_scripts/Environmental_Computing_ tutorials/data/Greenhouse.csv'
Greenhouse <- read.csv(file = Greenhouse_data, header = TRUE)

# response variable (Performance) to the left of the ~, the predictor variable (Treatment) to the right of the ~.
t.test(Performance ~ Treatment, data = Greenhouse, paired = TRUE)


############################################################
#                 Assumptions to check                     #
############################################################

# Normality
# For a paired t-test, it is assumed the the sample of differences is normally distributed. If these are highly skewed, 
# transformations may be used to achieve a distribution closer to normal.

# Independence
# The paired design takes into account that the two measures from each pair are not independent. It is still important, 
# however, that each pair of objects measured are independent from other pairs. If they are linked in any way (e.g., groups 
# of plants sharing a water tray) then more complex analytical design that account for additional factors may be required.


############################################################
#                Communicating the results                 #
############################################################

# Written
# As a minimum, the observed t statistic, the P-value and the number of degrees of freedom should be reported. 
# For example, you could write “The photosynthetic performance of plants was significantly greater in sunny environments 
# in contrast to shady environments (paired t-test: t = 18.81, df = 9, P < 0.001)”.

# Visual
# Box plots or column graphs with error bars are effective ways of communicating the variation in a single continuous 
# response variable versus a single categorical predictor variable.
boxplot(Performance~Treatment, data = Greenhouse, xlab = "Light environment", ylab = "Photosynthetic performance (FvFm)")

