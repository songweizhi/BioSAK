
################################# One_Sample_Wilcoxon_test #################################

# The one-sample Wilcoxon signed rank test is a non-parametric alternative to one-sample t-test when the data cannot be 
# assumed to be normally distributed. It’s used to determine whether the median of the sample is equal to a known standard 
# value (i.e. theoretical value).

# Note that, the data should be distributed symmetrically around the median. In other words, there should be roughly the 
# same number of values above and below the median.


# Here, we’ll use an example data set containing the weight of 10 mice. 
# We want to know, if the median weight of the mice differs from 25g?
set.seed(1234)
my_data <- data.frame(name = paste0(rep("M_", 10), 1:10), weight = round(rnorm(10, 20, 2), 1))
my_data


# test 
shapiro.test(my_data$weight)


# to be added


