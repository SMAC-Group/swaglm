rm(list=ls())
set.seed(1)
n <- 2000
p <- 20
Sigma <- diag(rep(1/p, p))

X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-10,5,6,19,70,rep(0,15))
z <- 1 + X%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))



fit1 = glm(y~X-1, family = binomial())
fit1
AIC(fit1)

fit2 = fastglm(x = X, y=as.numeric(y)-1, family=binomial())
test = summary(fit2)
test$coefficients[, 4]
summary(fit2)
fit2$coefficientsfit2$se
fit2$
?family
?fastglm
fit2
AIC(fit2)
fit2$aic



microbenchmark::microbenchmark(fit1 = glm(y~X-1, family = binomial()),
                               fit2 = fastglm(x = X, y=as.numeric(y)-1, family=binomial()))






