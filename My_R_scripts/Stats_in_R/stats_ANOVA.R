# ANOVA in R

y1 = c(18.2, 20.1, 17.6, 16.8, 18.8, 19.7, 19.1)
y2 = c(17.4, 18.7, 19.1, 16.4, 15.9, 18.4, 17.7)
y3 = c(15.2, 18.8, 17.7, 16.5, 15.9, 17.1, 16.7)

# combine all number list together
y = c(y1, y2, y3)

# generate label
n = rep(7, 3)
group = rep(1:3, n)

# Apply A Function Over An Array
tmp = tapply(y, group, stem)
stem(y)

tmpfn = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x))
tapply(y, group, tmpfn)

data = data.frame(y = y, group = factor(group))
data

# fit linear model
fit = lm(y ~ group, data)
fit
abline(fit)

plot(y ~ group, data)

# ANOVA
anova(fit)

?anova

