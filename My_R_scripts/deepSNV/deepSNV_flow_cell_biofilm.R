
# install GenomeInfoDbData first
source("https://bioconductor.org/biocLite.R")
biocLite("GenomeInfoDbData")

# update packages
update.packages('matrixStats')
update.packages('VGAM')

# install deepSNV (http://master.bioconductor.org/packages/devel/bioc/html/deepSNV.html), try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("deepSNV")

# load library
library(deepSNV)
?deepSNV

setwd('/Users/songweizhi/Desktop/deepSNV_bams')

test_bam = '2D9.bam'
control_bam = 'D2D0.bam'

# example 
regions = data.frame(chr="D2_p", start = 0, stop=1000000)

D2AD9 = deepSNV(test_bam, control_bam, regions=regions, q=25)
D2AD9

data(D2AD9)
show(D2AD9)
control(D2AD9)[100:110,]
test(D2AD9)[100:110,]
plot(D2AD9)
SNVs <- summary(D2AD9, sig.level=0.05, adjust.method="BH")
head(SNVs)
nrow(SNVs)
min(SNVs$freq.var)
sum(RF(test(D2AD9), total=T) > 0.01 & RF(test(D2AD9), total=T) < 0.95)
data(trueSNVs, package="deepSNV")
table(p.adjust(p.val(D2AD9), method="BH") < 0.05, trueSNVs)

######################################## Normalization ########################################

## Load data (unnormalized)
data(phiX, package="deepSNV")
plot(phiX, cex.min=.5)
## Normalize data
phiN <- normalize(phiX, round=TRUE) plot(phiN, cex.min=.5)
p.norm <- p.val(phiN)
## Normalize data
phiN <- normalize(phiX, round=TRUE) plot(phiN, cex.min=.5)
## Normalize data
phiN <- normalize(phiX, round=TRUE)
plot(phiN, cex.min=.5)
p.norm <- p.val(phiN)
n <- sum(!is.na(p.norm))
qqplot(p.norm, seq(1/n,1, length.out=n), log="xy", type="S", xlab="P-value", ylab="CDF")
p.val <- p.val(phiX)
points(sort(p.val[!is.na(p.val)]), seq(1/n,1, length.out=n), pch=16, col="grey", type="S", lty=2)
legend("topleft", c("raw data", "normalized data"), pch=16, col=c("grey", "black"), bty="n", lty= abline(0,1))
abline(0,1)

######################################## Overdispersion ########################################

data("RCC", package="deepSNV")
show(RCC)
plot(RCC, cex.min=.5)
RCC.bb = estimateDispersion(RCC, alternative="two.sided")
plot(RCC.bb, cex.min=.5)
RCC.bb@log.lik
RCC@log.lik
RCC.bb@log.lik - RCC@log.lik
log(4*nrow(test(RCC)))
summary(RCC, adjust.method="bonferroni")[,1:6]
tab <- summary(RCC.bb, adjust.method="bonferroni")[,1:6] tab
tab <- summary(RCC.bb, adjust.method="bonferroni")[,1:6]
tab
