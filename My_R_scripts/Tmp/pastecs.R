install.packages('pastecs')
library(pastecs)

identities = c(1,2,3,4,5,3,2,3,4,3,3,4)
abund(identities)
identities
plot(abund(identities))
data(bnr)
bnr
bnr.abd <- abund(bnr)
plot(bnr.abd, dpos=c(105, 100))
