# Introduction to mvabund
# http://environmentalcomputing.net/introduction-to-mvabund/

library(mvabund)


Herbivore_specialisation_data = '/Users/songweizhi/Desktop/COG_enrichment_analysis.csv'
Herbivore_specialisation_data = '/Users/songweizhi/Desktop/Herbivore_specialisation.csv'


Herbivores <- read.csv(file = Herbivore_specialisation_data, header = TRUE)
Herbivores


Herb_spp <- mvabund(Herbivores[,3:ncol(Herbivores)])
Herb_spp


par(mar=c(2,10,2,2)) # adjusts the margins
boxplot(Herbivores[,3:ncol(Herbivores)],horizontal = TRUE,las=2, main="Abundance")


meanvar.plot(Herb_spp)
plot(Herb_spp~Herbivores$Source, cex.axis=0.8, cex=0.8)
# use transformation="no" to look at the raw abundance data
# plot(Herb_spp~Herbivores$Habitat, cex.axis=0.8, cex=0.8, transformation="no")

# The model syntax below fits our response variable 
# (the mvabund object Herb_spp with the 100 counts of 7 species) 
# to the predictor variable Habitat (type of algae).

mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="poisson")
mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="negative.binomial")

# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="binomial")
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family=binomial(link="cloglog"))
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family=Gamma(link="log"))
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="binomial")
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="cloglog")
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="poisson")
# mod1 <- manyglm(Herb_spp ~ Herbivores$Source, family="gamma")

?manyglm
# If the model is a good fit, we should see a random scatter of points.
plot(mod1)



# We can test the multivariate hypothesis of whether species composition 
# varied across the habitats by using the anova function.
anova(mod1)
# We can see from this table that there is 
# a significant effect of Habitat (LRT = 625, P = 0.001),


# To examine this further, and see which herbivore species are more likely 
# to be found on which algal species, we can run univariate tests for each 
# species separately.
anova(mod1, p.uni="adjusted")

remove(list = ls())

