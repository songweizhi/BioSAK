
################################# Oneway ANOVA  #################################

# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

install.packages("devtools")
install.packages("ggpubr")
library("devtools")
library("ggpubr")
library("ggplot2")
library("dplyr")


################################# load data #################################

my_data <- PlantGrowth
my_data
str(my_data)
summary(my_data)
head(my_data)
View(my_data) # if you use RStudio this is a nice way to see the data in spreadsheet format
levels(my_data$group)


################################# Visualize data #################################

boxplot(my_data$weight)
boxplot(my_data$weight~my_data$group, col= rainbow(4))

# If the levels are not automatically in the correct order, re-order them as follow:
my_data$group <- ordered(my_data$group, levels = c("ctrl", "trt1", "trt2"))


ggboxplot(my_data, x = "group", y = "weight", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("ctrl", "trt1", "trt2"),
          ylab = "Weight", xlab = "Treatment")

ggline(my_data, x = "group", y = "weight", 
       add = c("mean_se", "jitter"), 
       order = c("ctrl", "trt1", "trt2"),
       ylab = "Weight", xlab = "Treatment")

# Box plot
boxplot(weight ~ group, data = my_data,
        xlab = "Treatment", ylab = "Weight",
        frame = FALSE, col = c("#00AFBB", "#E7B800", "#FC4E07"))

# plotmeans
install.packages("gplots")
library("gplots")
plotmeans(weight ~ group, data = my_data, frame = FALSE,
          xlab = "Treatment", ylab = "Weight",
          main="Mean Plot with 95% CI") 


################################# Compute one-way ANOVA test #################################

# Compute the analysis of variance
res.aov <- aov(weight ~ group, data = my_data)

# Summary of the analysis
summary(res.aov)

# As the p-value is less than the significance level 0.05, we can conclude that there are significant 
# differences between the groups highlighted with “*" in the model summary.


#################### Multiple pairwise-comparison between the means of groups ####################

# In one-way ANOVA test, a significant p-value indicates that some of the group means are different, 
# but we don’t know which pairs of groups are different.

# It’s possible to perform multiple pairwise-comparison, to determine if the mean difference between 
# specific pairs of group are statistically significant.


##### Tukey multiple pairwise-comparisons #####

# As the ANOVA test is significant, we can compute Tukey HSD (Tukey Honest Significant Differences) 
# for performing multiple pairwise-comparison between the means of groups.
# The function TukeyHD() takes the fitted ANOVA as an argument.
TukeyHSD(res.aov)

# It can be seen from the output, that only the difference between trt2 and trt1 is significant with 
# an adjusted p-value of 0.012.


##### Multiple comparisons using multcomp package #####

# It’s possible to use the function glht() [in multcomp package] to perform multiple comparison procedures for an ANOVA. 
# glht stands for general linear hypothesis tests. The simplified format is as follow:

# glht(model, lincft)
# 1) model: a fitted model, for example an object returned by aov().
# 2) lincft(): a specification of the linear hypotheses to be tested. Multiple comparisons in ANOVA models 
#    are specified by objects returned from the function mcp().


# perform multiple pairwise-comparisons for a one-way ANOVA:
library('multcomp')
summary(glht(res.aov, linfct = mcp(group = "Tukey")))


##### Pairewise t-test #####

# The function pairewise.t.test() can be also used to calculate pairwise comparisons between group levels 
# with corrections for multiple testing.
pairwise.t.test(my_data$weight, my_data$group, p.adjust.method = "BH")

# The result is a table of p-values for the pairwise comparisons. Here, the p-values have been adjusted 
# by the Benjamini-Hochberg method.


################################# Check ANOVA assumptions #################################

# The ANOVA test assumes that, the data are normally distributed and the variance across groups are homogeneous.
# We can check that with some diagnostic plots.


##### Check the homogeneity of variance #####

# The residuals versus fits plot can be used to check the homogeneity of variances
plot(res.aov, 1)
# In the plot, there is no evident relationships between residuals and fitted values (the mean of each groups), 
# We can assume the homogeneity of variances.
# Points 17, 15, 4 are detected as outliers, which can severely affect normality and homogeneity of variance. 
# It can be useful to remove outliers to meet the test assumptions.


# Levene’s test for testing homogeneity of variances.
# It’s also possible to use Bartlett’s test or Levene’s test to check the homogeneity of variances.
# We recommend Levene’s test, which is less sensitive to departures from normal distribution. 
library(car)
leveneTest(weight ~ group, data = my_data)
# From the output above we can see that the p-value is not less than the significance level of 0.05. 
# This means that there is no evidence to suggest that the variance across groups is statistically significantly different. 
# Therefore, we can assume the homogeneity of variances in the different treatment groups.


##### Relaxing the homogeneity of variance assumption #####

# The classical one-way ANOVA test requires an assumption of equal variances for all groups.
# How do we save our ANOVA test, in a situation where the homogeneity of variance assumption is violated?
# An alternative procedure (i.e.: Welch one-way test), that does not require that assumption have been 
# implemented in the function oneway.test().

# ANOVA test with no assumption of equal variances
oneway.test(weight ~ group, data = my_data)

# Pairwise t-tests with no assumption of equal variances
pairwise.t.test(my_data$weight, my_data$group, p.adjust.method = "BH", pool.sd = FALSE)


##### Check the normality assumption #####

# Normality plot of residuals
plot(res.aov, 2)

# In the plot below, the quantiles of the residuals are plotted against the quantiles of the normal distribution. 
# A 45-degree reference line is also plotted.

# The normal probability plot of residuals is used to check the assumption that the residuals are normally distributed. 
# It should approximately follow a straight line.

# As all the points fall approximately along this reference line, we can assume normality.


# Significance test with Shapiro-Wilk
# Extract the residuals
aov_residuals <- residuals(object = res.aov )

# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

