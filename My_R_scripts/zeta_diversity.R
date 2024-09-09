# zeta-diversity
library(zetadiv)

# Exploring zeta-diversity
# https://indradecastro.github.io/globalremoval/docs/20180525.zeta.html#read-data

# Measuring continuous compositional change using decline and decay in zeta diversity



utils::data(bird.spec.coarse)
xy.bird <- bird.spec.coarse[1:2]
data.spec.bird <- bird.spec.coarse[3:193]
dev.new(width = 12, height = 4)
zeta.bird <- Zeta.decline.mc(data.spec.bird, orders = 1:3, sam=100, plot = FALSE)
Plot.zeta.decline(zeta.bird)


zeta_demo_df <- read.csv(file = '/Users/songweizhi/Desktop/zeta_demo.txt', header = TRUE, sep = '\t')
value.demo = zeta_demo_df[2:ncol(zeta_demo_df)]
zeta.demo <- Zeta.decline.mc(value.demo, orders = 1:3, sam=100, plot = FALSE)
Plot.zeta.decline(zeta.demo)
