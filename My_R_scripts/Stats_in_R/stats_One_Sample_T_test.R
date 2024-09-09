
################################# One_Sample_T_test  #################################

# Link: http://www.sthda.com/english/wiki/one-sample-t-test-in-r


# Definition 
# one-sample t-test is used to compare the mean of one sample to a known standard (or theoretical/hypothetical) mean (μ).


# Generally, the theoretical mean comes from:
# 1) a previous experiment. For example, compare whether the mean weight of mice differs from 200 mg, a value determined 
#    in a previous study.
# 2) or from an experiment where you have control and treatment conditions. If you express your data as “percent of control”, 
#    you can test whether the average value of treatment condition differs significantly from 100.


# Note!!!
# One-sample t-test can be used only, when the data are normally distributed . This can be checked using Shapiro-Wilk test .


# Typical research questions are:
# 1) whether the mean (m) of the sample is equal to the theoretical mean (μ)?
# 2) whether the mean (m) of the sample is less than the theoretical mean (μ)?
# 3) whether the mean (m) of the sample is greater than the theoretical mean (μ)?

# In statistics, we can define the corresponding null hypothesis (H0) as follow:
# H0:m=μ
# H0:m≤μ
# H0:m≥μ

# The corresponding alternative hypotheses (Ha) are as follow:
# Ha:m≠μ (different)
# Ha:m>μ (greater)
# Ha:m<μ (less)

# If the p-value is inferior or equal to the significance level 0.05, we can reject the null hypothesis and accept 
# the alternative hypothesis. In other words, we conclude that the sample mean is significantly different from the 
# theoretical mean.


################################# Visualize data and compute one-sample t-test #################################

# install.packages("ggpubr")
# library('ggpubr')


# Here, we’ll use an example data set containing the weight of 10 mice, We want to know, if the average weight of the mice 
# differs from 25g?


# generate data
set.seed(1234)
my_data <- data.frame(name = paste0(rep("M_", 10), 1:10),weight = round(rnorm(10, 20, 2), 1))


# Check your data
head(my_data, 10) # print the first 10 rows
summary(my_data$weight) # Statistical summaries of weight
# Min.: the minimum value
# 1st Qu.: The first quartile. 25% of values are lower than this.
# Median: the median value. Half the values are lower; half are higher.
# 3rd Qu.: the third quartile. 75% of values are higher than this.
# Max.: the maximum value


# Visualize data with box plots
boxplot(my_data$weight)


################################# Preleminary test to check one-sample t-test assumptions #################################

# Is this a large sample? - No, because n < 30.
# Since the sample size is not large enough (less than 30, central limit theorem), we need to check whether the data 
# follow a normal distribution.


# check whether the data follow a normal distribution.
shapiro.test(my_data$weight)
# From the output, the p-value is greater than the significance level 0.05 implying that the distribution of the data 
# are not significantly different from normal distribtion. In other words, we can assume the normality.


# Visual inspection of the data normality using Q-Q plots (quantile-quantile plots). Q-Q plot draws the correlation 
# between a given sample and the normal distribution.
library("ggpubr")
ggqqplot(my_data$weight, ylab = "Weight", ggtheme = theme_minimal())

# From the normality plots, we conclude that the data may come from normal distributions.
# Note that, if the data are not normally distributed, it’s recommended to use the non parametric 
# one-sample Wilcoxon rank test.


################################# Compute one-sample t-test #################################

# One-sample t-test
res <- t.test(my_data$weight, mu = 25)
res 
# t is the t-test statistic value (t = -9.078),
# df is the degrees of freedom (df= 9),
# p-value is the significance level of the t-test (p-value = 7.95310^{-6}).
# conf.int is the confidence interval of the mean at 95% (conf.int = [17.8172, 20.6828]);
# sample estimates is he mean value of the sample (mean = 19.25).


#Access to the values returned by t.test() function
res$p.value # print the p-value
res$estimate # print the mean


# Test whether the mean weight of mice is less than 25g (one-tailed test)
res_2 = t.test(my_data$weight, mu = 25, alternative = "less")
res_2


# Test whether the mean weight of mice is greater than 25g (one-tailed test)
res_3 = t.test(my_data$weight, mu = 25, alternative = "greater")
res_3

