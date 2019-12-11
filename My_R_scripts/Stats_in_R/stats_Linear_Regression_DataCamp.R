
####################################### stat_Linear_Regression (DataCamp) #####################################

# https://www.datacamp.com/community/tutorials/linear-regression-R#coefficients

# Adding p values and R squared values to a plot using expression()
# http://lukemiller.org/index.php/2012/10/adding-p-values-and-r-squared-values-to-a-plot-using-expression/

#################################################### load data ################################################

library(readxl)
ageandheight <- read_excel("/Users/songweizhi/ageandheight.xls", sheet = "Hoja2") #Upload the data
ageandheight


########################################## Simple linear regression ###########################################

# Simple linear regression
lmHeight = lm(height~age, data = ageandheight) #Create the linear regression
summary(lmHeight) #Review the results


########################################## Multiple linear regression #########################################

# Multiple linear regression
lmHeight2 = lm(height~age + no_siblings, data = ageandheight) #Create a linear regression with two variables
summary(lmHeight2) #Review the results

# In general, for models that fit the data well, R² is near 1. Models that poorly fit the data have R² near 0.
# The adjusted R² is probably better to look at if you are adding more than one variable to the model


##################################################### Plot ####################################################

# plot the dots
plot(height~age, data = ageandheight, pch = 16, col = "blue") #Plot the results


# add the regression line
abline(lmHeight)


# get the value of r2 and p
regression_summary = summary(lmHeight)
r2 = regression_summary$adj.r.squared
p_value = regression_summary$coefficients[2,4]


# add r2 and p to the plot
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), list(MYOTHERVALUE = format(p_value, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n')


# plot residuals
# Ideally, when you plot the residuals, they should look random. 
# Otherwise means that maybe there is a hidden pattern that the linear model is not considering.
plot(lmHeight$residuals, pch = 16, col = "red")


