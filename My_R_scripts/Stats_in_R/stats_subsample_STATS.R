library(vegan)
library(dunn.test)


Sample = c(30,40,10,5,15)


df = NULL
for(i in 1:100){
  df = rbind(df, rrarefy(x = Sample, sample = 70))
}

df

boxplot(df)


G = c(rep("A", 100),rep("B", 100),rep("C", 100),rep("D", 100),rep("E", 100))
G
X = c(df[,1], df[,2], df[,3], df[,4], df[,5])
X

AOV = aov(X ~ G)
summary(AOV)
pairwise.t.test(x = X, g = G, p.adjust.method = "BH")

dunn.test(x = X, g = G)


