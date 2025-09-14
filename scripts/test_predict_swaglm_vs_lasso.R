rm(list=ls())



# ----------------------------------- Simulation comparing SWAG and Lasso for prediction


# Load required package
library(glmnet)
library(swaglm)
library(pROC)

# Example data: train/test split
n <- 2000
p <- 100
# create design matrix and vector of coefficients
Sigma <- diag(rep(1/p, p))
seed= 12345
set.seed(seed)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-15,-10,5,10,15, rep(0,p-5))

# --------------------- generate from logistic regression with an intercept of one


B = 500
names_mat = c("seed",
              "lambda","train_accuracy_lasso", "test_accuracy_lasso", "train_auc_lasso", "test_auc_lasso","time_lasso",
              "train_accuracy_swaglm", "test_accuracy_swaglm", "train_auc_swaglm", "test_auc_swaglm", "time_swaglm")
mat_res = matrix(NA, ncol=length(names_mat), nrow=B)
colnames(mat_res) = names_mat


for(b in seq(B)){
  # b=1
  seed = 12345 + b
  z <- 1 + X%*%beta
  pr <- 1/(1 + exp(-z))
  set.seed(seed)
  y <- as.factor(rbinom(n, 1, pr))
  y = as.numeric(y)-1
  
  # Split into train/test
  prop_train = 0.7
  train_idx <- sample(seq_len(n), size = round(prop_train * n))
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test  <- X[-train_idx, ]
  y_test  <- y[-train_idx]
  
  # ---------------------------- Fit LASSO logistic regression with cross-validation
  
  t1 = Sys.time()
  cvfit <- cv.glmnet(
    x = X_train, 
    y = y_train,
    family = "binomial",
    alpha = 1,          # LASSO penalty
    # type.measure = "class" # can also use "deviance" or "auc"
  )
  
  # Best lambda
  lambda_opt <- cvfit$lambda.min
  
  # Predict on train and test
  prob_train_lasso <- predict(cvfit, newx = X_train, s = "lambda.min", type = "response")
  prob_test_lasso  <- predict(cvfit, newx = X_test,  s = "lambda.min", type = "response")
  
  # Class predictions
  pred_train_lasso <- ifelse(prob_train_lasso > 0.5, 1, 0)
  pred_test_lasso  <- ifelse(prob_test_lasso > 0.5, 1, 0)
  
  t2 = Sys.time()
  
  time_lasso = difftime(t2,t1,units = "secs")
  
  # accuracy
  acc_train_lasso <- mean(pred_train_lasso == y_train)
  acc_test_lasso  <- mean(pred_test_lasso == y_test)
  

  
  # AUC 
  auc_train_lasso <- suppressMessages(auc(y_train, as.numeric(pred_train_lasso)))
  auc_test_lasso <- suppressMessages(auc(y_test, as.numeric(pred_test_lasso)))
  
  # save in mat
  mat_res[b,"lambda" ] = lambda_opt
  mat_res[b, "train_accuracy_lasso"] = acc_train_lasso
  mat_res[b,"test_accuracy_lasso" ] = acc_test_lasso
  mat_res[b,"train_auc_lasso" ] =   auc_train_lasso
  mat_res[b,"test_auc_lasso" ] =   auc_test_lasso
  mat_res[b,"time_lasso" ] =   as.numeric(time_lasso)
  
  # ------------ swaglm
  
  
  t3 = Sys.time()
  alpha=0.15
  pmax=10
  fit_swaglm = swaglm(X_train, y_train, p_max = pmax, family = binomial(),
                      alpha =alpha, verbose = FALSE)
  y_pred_train_swaglm = predict(fit_swaglm, X_train)  
  y_pred_test_swaglm =  predict(fit_swaglm, X_test)
  
  #
  y_pred_train_swaglm_majority_class = ifelse(rowMeans(y_pred_train_swaglm$mat_reponse_prediction)>.5, 1, 0)
  y_pred_test_swaglm_majority_class = ifelse(rowMeans(y_pred_test_swaglm$mat_reponse_prediction)>.5, 1, 0)
  
  y_pred_train_swaglm_mean_pred = rowMeans(y_pred_train_swaglm$mat_reponse_prediction)
  y_pred_test_swaglm_mean_pred = rowMeans(y_pred_test_swaglm$mat_reponse_prediction)
  
  
  t4 = Sys.time()
  time_swaglm = difftime(t4,t3,units = "secs")
  
  
  # accuracy
  acc_train_swaglm <- mean(y_pred_train_swaglm_majority_class == y_train)
  acc_test_swaglm  <- mean(y_pred_test_swaglm_majority_class == y_test)
  

  
  # AUC 
  auc_train_swaglm <- suppressMessages(auc(y_train, as.numeric(y_pred_train_swaglm_mean_pred)))
  auc_test_swaglm  <- suppressMessages(auc(y_test, as.numeric(y_pred_test_swaglm_mean_pred)))
  
  mat_res[b, "train_accuracy_swaglm"] = acc_train_swaglm
  mat_res[b,"test_accuracy_swaglm" ] = acc_test_swaglm
  mat_res[b,"train_auc_swaglm" ] =   auc_train_swaglm
  mat_res[b,"test_auc_swaglm" ] =   auc_test_swaglm
  mat_res[b,"time_swaglm" ] =   time_swaglm
  
  cat(paste0(b, " \n"))
  }



# Store results



# print(results)