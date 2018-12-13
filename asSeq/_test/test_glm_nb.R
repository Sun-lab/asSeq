library(MASS)
library(asSeq)

x1 = rep(0:2, times=c(100, 200, 100))
x2 = sample(x1)

y  = rep(NA, 400)
mu = 2 + 0.1*x1 + 0.1*x2

for(i in 1:400){
  y[i] = rnegbin(1, exp(mu[i]), theta=2)
}

# quartz(width=8, height=4)
par(mfrow=c(1,2))
boxplot(y~x1)
boxplot(y~x1+x2)

g1 = glm.nb(y~x1 + x2)
summary(g1)

X  = cbind(x1, x2)
g2 = glmNB(y, X, trace=1)

summary(g1$resid - g2$resid)

g1$theta
1/g2$phi

logLik(g1)
g2$loglik
