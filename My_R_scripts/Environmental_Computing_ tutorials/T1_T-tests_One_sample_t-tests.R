# One sample t-tests
# http://environmentalcomputing.net/one-sample-t-test/

pollen <- c(94,135,78,98,137,114,114,101,112,121)
t.test(pollen, mu = 125)


############################################################
#                 Assumptions to check                     #
############################################################

# Normality
# The t distribution describes paramaters sampled from a normal population, so assumes 
# that the data are normally distributed. Note however that t-tests are reasonably robust 
# to violations of normality (although watch out for the influence of outliers).

# Independence
# The observations should have been sampled randomly from a defined population so that 
# sample mean is an unbiased estimate of the population mean. If individual replicates 
# are linked in any way, then the assumption of independence will be violated.


############################################################
#                Communicating the results                 #
############################################################

# Written 
# As a minimum, the observed t statistic, the p value and the number of degrees of freedom 
# should be reported. 
# For example, you could write “The mean pollen count from the footprints (109 ) was 
# significantly lower than expected if it was derived from the nearby forest with an 
# average count of 125 (t = 3.07, df = 9, P = 0.01)”.

# Visual
# Box plots or frequency histograms can be used to visualise the variation in a single variable. 
# In this example, you might use a line or arrow to indicate the single value (125) that you were 
# comparing the sample to.

hist(pollen,main=NULL)
abline(v=125,col="red")

boxplot(pollen,xlab = "Pollen count",horizontal =TRUE)
abline(v=125,col="red")

