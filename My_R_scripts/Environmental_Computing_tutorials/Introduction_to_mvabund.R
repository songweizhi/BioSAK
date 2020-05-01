# Introduction to mvabund
# http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)


Herbivore_specialisation_data = '/Users/songweizhi/PycharmProjects/BioSAK/My_R_scripts/Environmental_Computing_ tutorials/data/Herbivore_specialisation.csv'


Herbivores <- read.csv(file = Herbivore_specialisation_data, header = TRUE)
Herbivores


Herb_spp <- mvabund(Herbivores[,5:11])
Herb_spp


par(mar=c(2,10,2,2)) # adjusts the margins
boxplot(Herbivores[,5:11],horizontal = TRUE,las=2, main="Abundance")


meanvar.plot(Herb_spp)


plot(Herb_spp~Herbivores$Habitat, cex.axis=0.8, cex=0.8)
# use transformation="no" to look at the raw abundance data
# plot(Herb_spp~Herbivores$Habitat, cex.axis=0.8, cex=0.8, transformation="no")

# The model syntax below fits our response variable 
# (the mvabund object Herb_spp with the 100 counts of 7 species) 
# to the predictor variable Habitat (type of algae).
mod1 <- manyglm(Herb_spp ~ Herbivores$Habitat, family="poisson")

# If the model is a good fit, we should see a random scatter of points.
plot(mod1)


# We can use the family argument to choose a distribution 
# which is better suited to our data. 
mod2 <- manyglm(Herb_spp ~ Herbivores$Habitat, family="negative_binomial")


# This residual plot is much better, there is now no discernible fan shape 
# and we will use this model for all further analysis.
plot(mod2)


# We can test the multivariate hypothesis of whether species composition 
# varied across the habitats by using the anova function.
anova(mod2)
# We can see from this table that there is 
# a significant effect of Habitat (LRT = 625, P = 0.001),


# To examine this further, and see which herbivore species are more likely 
# to be found on which algal species, we can run univariate tests for each 
# species separately.
anova(mod2, p.uni="adjusted")


# we can test more complex models. For example, to fit a model with both habitat and day or night
mod3 <- manyglm(Herb_spp ~ Herbivores$Habitat*Herbivores$DayNight, family="negative_binomial")
anova(mod3)
# You can see that the species composition of herbivores varies with habitat, but not between day and night.



