
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

# example 
regions <- data.frame(chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 3120, stop=3140)
HIVmix100 <- deepSNV(test = system.file("extdata", "test.bam", package="deepSNV"), control = system.file("extdata", "control.bam", package="deepSNV"), regions=regions, q=10)

HIVmix200 <- deepSNV(test = "/Users/songweizhi/Desktop/deepSNV_bams/test.bam", 
                   control = "/Users/songweizhi/Desktop/deepSNV_bams/control.bam", 
                   regions=regions, 
                   q=10)
HIVmix100

HIVmix200

#HIVmix2 <- deepSNV(test = '/Users/songweizhi/Desktop/deepSNV_bams/test.bam', control = '/Users/songweizhi/Desktop/deepSNV_bams/control.bam', regions=regions, q=10)
HIVmix3
data(HIVmix3)
show(HIVmix3)
show(HIVmix3)
control(HIVmix2)[100:110,]
test(HIVmix2)[100:110,]
plot(HIVmix2)



data(HIVmix)
show(HIVmix)
show(HIVmix)
control(HIVmix)[100:110,]
test(HIVmix)[100:110,]
plot(HIVmix)
SNVs <- summary(HIVmix, sig.level=0.05, adjust.method="BH")
head(SNVs)
nrow(SNVs)
min(SNVs$freq.var)
sum(RF(test(HIVmix), total=T) > 0.01 & RF(test(HIVmix), total=T) < 0.95)
data(trueSNVs, package="deepSNV")
table(p.adjust(p.val(HIVmix), method="BH") < 0.05, trueSNVs)

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
