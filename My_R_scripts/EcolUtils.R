
library(devtools)
devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)


?spec.gen
#spec.gen(comm.tab, niche.width.method = "levins", perm.method = "quasiswap", n = 1000, probs = c(0.025, 0.975))


?as.numeric
??niche.width
??occurrence

# niche.width.method=="occurrence"
levin.index.real<-occurrence(comm.tab)

# niche.width.method=="levins"
levin.index.real<-as.numeric(niche.width(comm.tab,method="levins"))

# niche.width.method=="shannon"
levin.index.real<-as.numeric(niche.width(comm.tab,method="shannon"))

# nee_niche_breadth
levin.index.real<-nee_niche_breadth(comm.tab)


