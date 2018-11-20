
library(asSeq)

# -------------------------------------------------------------
# simulate data
# -------------------------------------------------------------

x1 = rep(0:2, times=c(100, 200, 100))
x2 = sample(x1)
x3 = x1 + rnorm(400, 0, 0.5)

y  = rep(NA, 400)
mu = 2 + 0.2*x1 + 0.1*x2 + 0.1*x3

for(i in 1:400){
  y[i] = rpois(1, exp(mu[i]))
}

boxplot(y~x1+x2)

# -------------------------------------------------------------
# use R/glm
# -------------------------------------------------------------

g1 = glm(y~x1 + x2 + x3, family=poisson)
summary(g1)

# -------------------------------------------------------------
# use glmFit
# -------------------------------------------------------------

X  = cbind(x1, x2, x3)
g2 = glmFit(family="poisson", link="log", y, X)

# -------------------------------------------------------------
# compare residuals
# -------------------------------------------------------------

if( max(abs(g1$resid - g2$resid))> 1e-6 ){
  stop("residuals are different\n")
}

summary(g1$resid - g2$resid)

# -------------------------------------------------------------
# test glmTest
# -------------------------------------------------------------

X  = cbind(x1, x2)
gm = glmFit(family="poisson", link="log", y, X)

Z  = matrix(x3, ncol=1)
glmTest(gm, Z)

g1 = glm(y~x1 + x2, family=poisson)
g2 = glm(y~x1 + x2 + x3, family=poisson)
anova(g1, g2, test="Chisq")
