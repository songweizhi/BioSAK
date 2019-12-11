install.packages("ggplot2")
library(ggplot2)

d
l_l = d
jpeg(filename = '/Users/weizhisong/Desktop/L_L.jpg')
hist(d, axes = T, xlab = 'Identity', xlim = range(60:85),ylab = 'Target Frequency', probability = T, main = 'Group: L_L', border = 'blue')
dev.off()
lines(density(l_l,width = 5), lwd = 2)
dev.off()

e = c(1,2,3,4,5,4,3,2,3,3,4,3,4,3,4,5,6)
hist(e)
dev.off()
lines(density(e,width = 5), lwd = 2)
dev.off()



