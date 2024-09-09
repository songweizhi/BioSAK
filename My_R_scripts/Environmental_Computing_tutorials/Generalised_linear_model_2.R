# Generalised linear models 2
# http://environmentalcomputing.net/generalised-linear-models-2/

library(mvabund)

Revegetation_data = '/Users/songweizhi/PycharmProjects/BioSAK/My_R_scripts/Environmental_Computing_tutorials/data/Revegetation.csv'
Reveg <- read.csv(Revegetation_data, header = T)
Reveg


hist(Reveg$Soleolifera)


ft.sol.pois <- manyglm(Soleolifera~Treatment, family="poisson", data=Reveg)
plot(ft.sol.pois)


ft.sol.nb <- manyglm(Soleolifera~Treatment,family="negative binomial",data=Reveg)
plot(ft.sol.nb)


anova(ft.sol.nb)
summary(ft.sol.nb)
ft.sol.nb$coefficients

boxplot(Soleolifera~Treatment,ylab = "Count", xlab = "Treatment", data=Reveg)



