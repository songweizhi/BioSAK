faithful
duration = faithful$eruptions
duration
quantile(duration, c(.05, .32, .57, .98)) 

plot(density(duration))
barplot(duration)

help(quantile)
write.csv(duration, file = '/Users/weizhisong/Desktop/duration.csv')
