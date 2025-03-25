library(speedglm)
library(microbenchmark)
library(ggplot2)
library(fastglm)
# install.packages("speedglm")

set.seed(123)
n.obs  <- 10000
n.vars <- 100
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
Sigma <- 0.99 ^ abs(outer(1:n.vars, 1:n.vars, FUN = "-"))
x <- MASS::mvrnorm(n.obs, mu = runif(n.vars, min = -1), Sigma = Sigma)

y <- 1 * ( drop(x[,1:25] %*% runif(25, min = -0.1, max = 0.10)) > rnorm(n.obs))

ct <- microbenchmark(
  glm.fit = {gl1 <- glm.fit(x, y, family = binomial())},
  speedglm.eigen  = {sg1 <- speedglm.wfit(y, x, intercept = FALSE,
                                          family = binomial())},
  speedglm.chol   = {sg2 <- speedglm.wfit(y, x, intercept = FALSE,
                                          family = binomial(), method = "Chol")},
  speedglm.qr     = {sg3 <- speedglm.wfit(y, x, intercept = FALSE,
                                          family = binomial(), method = "qr")},
  fastglm.qr.cpiv = {gf1 <- fastglm(x, y, family = binomial())},
  fastglm.qr      = {gf2 <- fastglm(x, y, family = binomial(), method = 1)},
  fastglm.LLT     = {gf3 <- fastglm(x, y, family = binomial(), method = 2)},
  fastglm.LDLT    = {gf4 <- fastglm(x, y, family = binomial(), method = 3)},
  fastglm.qr.fpiv = {gf5 <- fastglm(x, y, family = binomial(), method = 4)},
  times = 25L
)

autoplot(ct, log = FALSE) + stat_summary(fun.y = median, geom = 'point', size = 2)
