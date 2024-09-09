
################################# stat_Linear_Regression #################################

# http://r-statistics.co/Linear-Regression.html


######################################## load data #######################################

cars


###################################### Scatter Plot ######################################

scatter.smooth(x=cars$speed, y=cars$dist, main="Dist ~ Speed")  # scatterplot


############################## BoxPlot – Check for outliers ##############################

# divide graph area in 2 columns
par(mfrow=c(1, 2))  

# box plot for 'speed'
boxplot(cars$speed, main="Speed", sub=paste("Outlier rows: ", boxplot.stats(cars$speed)$out))  

# box plot for 'distance'
boxplot(cars$dist, main="Distance", sub=paste("Outlier rows: ", boxplot.stats(cars$dist)$out))  


###################################### Density plot ######################################

# To check if the response variable is close to normality

# load library
library(e1071)

# divide graph area in 2 columns
par(mfrow=c(1, 2))  

# density plot for 'speed'
plot(density(cars$speed), main="Density Plot: Speed", ylab="Frequency", sub=paste("Skewness:", round(e1071::skewness(cars$speed), 2)))  
polygon(density(cars$speed), col="red")

# density plot for 'dist'
plot(density(cars$dist), main="Density Plot: Distance", ylab="Frequency", sub=paste("Skewness:", round(e1071::skewness(cars$dist), 2)))  
polygon(density(cars$dist), col="red")


####################################### Correlation ######################################

# calculate correlation between speed and distance 
cor(cars$speed, cars$dist)  


################################### Build Linear Model ###################################

# build linear regression model on full data
linearMod <- lm(dist ~ speed, data=cars)  
print(linearMod)


############################## Linear Regression Diagnostics #############################

# print the summary statistics for obtained linear model
modelSummary <- summary(linearMod)
modelSummary

# model p-Value (bottom last line: p-value: 1.49e-12)
# p-Value of individual predictor variables (extreme right column under ‘Coefficients’)

# The p-Values are very important because, We can consider a linear model to be statistically 
# significant only when both these p-Values are less that the pre-determined statistical 
# significance level, which is ideally 0.05.

# get individual values

# get model coefficients
modelCoeffs <- modelSummary$coefficients 
modelCoeffs

























