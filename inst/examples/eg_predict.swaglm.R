set.seed(12345)
n <- 2000
p <- 100
# create design matrix and vector of coefficients
Sigma <- diag(rep(1/p, p))
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(10, rep(0,p-1))
sigma2=2
y <- 1 + X%*%beta + rnorm(n, mean = 0, sd = sqrt(sigma2))

# subset data
n_data_train = n-5
X_sub = X[1:n_data_train, ]
y_sub = y[1:n_data_train]

# plot train data
plot(X_sub[,1], y_sub)
abline(a=1, b=beta[1])

# run swag
swaglm_obj = swaglm(X = X_sub, y = y_sub, p_max = 15, family = gaussian(),
                    method = 0, alpha = .15, verbose = TRUE)

# compute prediction
X_to_predict = X[(n_data_train+1):(dim(X)[1]), ]
y_pred = predict(swaglm_obj, newdata = X_to_predict)
n_selected_model = dim(y_pred$mat_eta_prediction)[2]

# in that case mat_eta_prediction and mat_reponse_prediction are the same
all.equal(y_pred$mat_eta_prediction, y_pred$mat_reponse_prediction)

# compute average prediction (accross selected models)
y_pred_mean = apply(y_pred$mat_reponse_prediction, MARGIN = 1, FUN = mean)

# plot
y_to_predict = y[(n_data_train+1):(dim(X)[1])]
plot(X_to_predict[,1], y_to_predict, col="red", ylim=c(-4,4), ylab="y", xlab="X")
# add all prediction
for(i in seq(dim(X_to_predict)[1])){
  # i=1
  points(x = rep(X_to_predict[i,1],n_selected_model ),
         y = y_pred$mat_reponse_prediction[i,])
}
points(x =X_to_predict[,1], y = y_pred_mean, col="orange")
abline(a=1, b=beta[1])
legend(
  "topleft",
  legend = c(
    "True values (test set)",
    "Predictions (all selected models)",
    "Average prediction",
    "True regression line"
  ),
  col = c("red", "black", "orange", "black"),
  pch = c(1, 1, 1, NA),
  lty = c(NA, NA, NA, 1),
  cex = 0.8,
  bty="n"
)



# ----------------------------------- logistic regression

set.seed(12345)
n <- 2000
p <- 100
# create design matrix and vector of coefficients
Sigma <- diag(rep(1/p, p))
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-15,-10,5,10,15, rep(0,p-5))

# --------------------- generate from logistic regression with an intercept of one
z <- 1 + X%*%beta
pr <- 1/(1 + exp(-z))
y <- as.factor(rbinom(n, 1, pr))
y = as.numeric(y)-1

# subset data
n_data_train = n-200
X_sub = X[1:n_data_train, ]
y_sub = y[1:n_data_train]

swaglm_obj = swaglm::swaglm(X=X_sub, y = y_sub, p_max = 20, family = stats::binomial(),
                            alpha = .15, verbose = TRUE, seed = 123)

# compute prediction
X_to_predict = X[(n_data_train+1):(dim(X)[1]), ]
y_pred = predict(swaglm_obj, newdata = X_to_predict)
y_pred_majority_class <- ifelse(rowMeans(y_pred$mat_reponse_prediction) >= 0.5, 1, 0)
y_to_predict = y[(n_data_train+1):(dim(X)[1])]

# tabulate
table(True = y_to_predict, Predicted = y_pred_majority_class)

