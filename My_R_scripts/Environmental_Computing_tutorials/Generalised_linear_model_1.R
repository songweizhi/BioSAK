library(mvabund)

# Generalised linear models 1
# http://environmentalcomputing.net/generalised-linear-models-1/

Crabs_data = '/Users/songweizhi/PycharmProjects/BioSAK/My_R_scripts/Environmental_Computing_tutorials/data/Crabs.csv'

Crab_PA <- read.csv(Crabs_data, header = T)
Crab_PA

# To test whether the probability of crab presence changes with time (a factor) 
# and distance (a continuous variable), we fit the following model. 
# The response variable (presence/absence of crabs) is binomial, so we use family=binomial.
ft.crab <- glm(CrabPres ~ Time*Dist, family=binomial, data=Crab_PA)

# Assumptions of GLM
# Assumption 1 : The observed y are independent, conditional on some predictors x.
# Assumption 2 : The response y come from a known distribution with a known mean-variance relationship.
# To check Assumption 2, we look at a plot of residuals, and try to see if there is a fan shape.
plot(ft.crab,which=1)
# Unfortunately the glm plot function gives us a very odd looking plot due to the discreteness of the data.

# For a more useful plot we can instead fit the model using the manyglm function in the mvabund package.
ft.crab.many <- manyglm(CrabPres ~ Time*Dist, family="binomial", data=Crab_PA)
plot(ft.crab.many)
# there doesn’t seem to be a fan shape, so we can conclude the mean-variance assumption the model made was reasonable for our data.
# The residuals in this plot have a random component. If you see a pattern it’s best to repeat the plot a few times to see if the pattern is real.

# Assumption 3 : There is a straight line relationship between a known function g of the mean of y and the predictors x
# To check this assumption, we check the residual plot above for non-linearity, 
# or a U-shape. In our case there is no evidence of non-linearity. 
# If the residuals seem to go down then up, or up then down, we may need to 
# add a polynomial function of the predictors using the poly function.


# Interpreting the results
anova(ft.crab,test="Chisq")

anova(ft.crab.many)

summary(ft.crab)









