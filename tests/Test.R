# Tests #
library(convReg)
##### Binom + Gauss #####
n = 200
set.seed(123)
x1 = rnorm(n,3,0.5)
x2 = runif(n,-2,1)
k= rbinom(n,prob=exp(-exp(x2)),size=2)
e=x1+rnorm(n=n,mean=0,sd=0.25)
y= data.frame(obs=(k + e) , f1 = x1, f2 = x2)
par(mfrow = c(1,1), mar = c(5,4,4,2))
hist(y$obs, breaks = seq(min(y$obs)-1, max(y$obs)+1, length.out = n/20), xlim = c(min(y$obs)-1, max(y$obs)+1))
dist1 = "Binom";dist2 = "Gauss"
fixed=list(name = "sigma 1: (Intercept)", value = 2)
res.reg.em =convreg( ~obs,
           formula.mu1 =~ f2,
           formula.mu2 =~ f1,
           fixed = fixed,
           data=y,dist1 = dist1,
           method = "em",
           scale=T, scaleInit = max(y$obs)/2)

bS = BICselect(df = y, formula.resp = ~obs, idx.pred = 2:3, dist1 = "Binom", dist2 = "Gauss")

td = testdist(y$obs, dist1 = "Binom", dist2 = "Gauss",fixed=list(name = "sigma 1: (Intercept)", value = 3),method = "em",
              scale=T, scaleInit = max(y$obs)/3)
plotTestdist(td)

td = testdist(y$obs, dist1 = "Pois", dist2 = "Lnorm")
plotTestdist(td)

summary(res.reg.em)
plot(res.reg.em)
regplot(res.reg.em)
distplot(res.reg.em)
cdfplot(res.reg.em)


####

x1 = rnorm(n,3,0.5)
x2 = runif(n,-2,1)
k= rpois(n,exp(-exp(x2)))
e=x1+rlnorm(n=n,mean=0,sd=0.25)
y= data.frame(obs=(k + e) , f1 = x1, f2 = x2)
par(mfrow = c(1,1), mar = c(5,4,4,2))
hist(y$obs, breaks = seq(min(y$obs)-1, max(y$obs)+1, length.out = n/5), xlim = c(min(y$obs)-1, max(y$obs)+1))
dist1 = "Pois";dist2 = "Lnorm"
res.reg.em =convreg( ~obs,
                     formula.mu1 =~ f2,
                     formula.mu2 =~ f1,
                     data=y,
                     dist1 = dist1,
                     dist2 = dist2,
                     method = "mle")

bS = BICselect(df = y, formula.resp = ~obs, idx.pred = 2:3, dist1 = dist1, dist2 = dist2)

td = testdist(y$obs, dist1 = dist1, dist2 = dist2)
plotTestdist(td)

td = testdist(y$obs, dist1 = dist1, dist2 = "Gauss")
plotTestdist(td)

summary(res.reg.em)
plot(res.reg.em)
regplot(res.reg.em)
distplot(res.reg.em)
cdfplot(res.reg.em)

