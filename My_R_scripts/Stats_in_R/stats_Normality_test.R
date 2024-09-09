
################################# Normality Test  #################################

# Many of statistical tests including correlation, regression, t-test, and analysis of variance (ANOVA) assume some 
# certain characteristics about the data. They require the data to follow a normal distribution or Gaussian distribution. 
# These tests are called parametric tests, because their validity depends on the distribution of the data.

# Before using a parametric test, we should perform some preleminary tests to make sure that the test assumptions are met. 
# In the situations where the assumptions are violated, non-paramatric tests are recommended.

# Here, we’ll describe how to check the normality of the data by visual inspection and by significance tests.


################################# Load packages and Import data #################################

# Install required R packages
#install.packages("dplyr")
#install.packages("ggpubr")
library("dplyr")
library("ggpubr")


# Import data
my_data <- ToothGrowth
my_data


# print a random sample of 10 rows using the sample_n() function
set.seed(1234)
dplyr::sample_n(my_data, 10)


# Note!!!
# Case of large sample sizes:
# If the sample size is large enough (n > 30), we can ignore the distribution of the data and use parametric tests.
# The central limit theorem tells us that no matter what distribution things have, the sampling distribution tends 
# to be normal if the sample is large enough (n > 30).


# Normality can be checked by visual inspection [normal plots (histogram), Q-Q plot (quantile-quantile plot)] or by 
# significance tests].


################################# Normality test (Visual methods) #################################

# Visual methods: Density plot
library("ggpubr")
ggdensity(my_data$len, main = "Density plot of tooth length", xlab = "Tooth length")


# Visual methods: Q-Q plot
library("ggpubr")
ggqqplot(my_data$len)
# As all the points fall approximately along this reference line, we can assume normality.


################################# Normality test (Significance test) #################################

# Visual inspection is usually unreliable. It’s possible to use a significance test comparing the sample distribution 
# to a normal one in order to ascertain whether data show or not a serious deviation from normality.

# Shapiro-Wilk’s method is widely recommended for normality test and it provides better power than K-S. It is based on the 
# correlation between the data and the corresponding normal scores.

# The null hypothesis of these tests is that “sample distribution is normal”. If the test is significant, the distribution 
# is non-normal.

# Note that, normality test is sensitive to sample size. Small samples most often pass normality tests. Therefore, it’s 
# important to combine visual inspection and significance test in order to take the right decision.


# perform Shapiro-Wilk test
shapiro.test(my_data$len)
#  The p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. 
# In other words, we can assume the normality.

