
library(mvabund)

dat  <- read.delim("/Users/songweizhi/Desktop/test/Pathways_PA_energy.txt", row.names = 1)
fact <- read.delim("/Users/songweizhi/Desktop/test/genome_cate.txt", row.names = 1)

dat.mv <- mvabund(dat)
dat.nb = manyglm(dat.mv ~ Group, family="binomial", data=fact)

# Check the model fitting plotting the residuals

plot(dat.nb,n.vars=15)


# run the statistical test
dat.aov = anova(dat.nb, nBoot=1, p.uni="unadjusted") # change to 1000

dat.aov$uni.p

# multiple comparison
# here fit the pathway of interest as GLM (one allowed each time so if you have multiple to do, try for loop)
# for example
PWY1YI0.3 <- dat$PWY1YI0.3

# make mvabund subject
PWY1YI0.3 <- mvabund(PWY1YI0.3)

ft.1 <- manyglm(PWY1YI0.3 ~ Group, family="binomial", data=fact)

plot(ft.1,n.vars=15)


anova(ft.1 , pairwise.comp = fact$Group, nBoot = 1000)

# results

Time elapsed: 0 hr 0 min 1 sec
Analysis of Deviance Table

Model: PWY1YI0.3 ~ Group

Multivariate test:
  Res.Df Df.diff   Dev Pr(>Dev)
(Intercept)   1226                       
Group         1224       2 2.267    0.307

Pairwise comparison results: 
  Observed statistic Free Stepdown Adjusted P-Value
BothAbsent vs BothFound              1.794                          0.457
BothFound vs OnlyNif                 1.340                          0.457
BothAbsent vs OnlyNif                0.009                          0.910

Arguments: P-value calculated using 1000 iterations via PIT-trap resampling.



