rm(list=ls())

# ----------------------------------- Simulation comparing SWAG and Lasso for prediction
library(glmnet)
library(swaglm)
library(pROC)
library(MASS)        # for mvrnorm
library(doParallel)  # for parallel execution
library(foreach)
library(beepr)

# Example data: train/test split
n <- 1000
p <- 100
Sigma <- diag(rep(1/p, p))
seed= 12345
set.seed(seed)
X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
beta = c(-15,-10,5,10,15, rep(0,p-5))

# --------------------- simulation setup
B = 100
names_mat = c("seed",
              "lambda","train_accuracy_lasso", "test_accuracy_lasso", "train_auc_lasso", "test_auc_lasso","time_lasso",
              "train_accuracy_swaglm", "test_accuracy_swaglm", "train_auc_swaglm", "test_auc_swaglm", "time_swaglm")

# Set up cluster
ncores <- parallel::detectCores() - 1
ncores
cl <- makeCluster(ncores)
registerDoParallel(cl)

# Run parallel Monte Carlo
mat_res <- foreach(b = seq(B), .combine = rbind,
                   .packages = c("glmnet","swaglm","pROC","MASS")) %dopar% {
                     
                     # seed for reproducibility
                     seed = 12345 + b
                     z <- 1 + X%*%beta
                     pr <- 1/(1 + exp(-z))
                     set.seed(seed)
                     y <- rbinom(n, 1, pr)
                     
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
                       alpha = 1
                     )
                     
                     lambda_opt <- cvfit$lambda.min
                     
                     prob_train_lasso <- predict(cvfit, newx = X_train, s = "lambda.min", type = "response")
                     prob_test_lasso  <- predict(cvfit, newx = X_test,  s = "lambda.min", type = "response")
                     
                     pred_train_lasso <- ifelse(prob_train_lasso > 0.5, 1, 0)
                     pred_test_lasso  <- ifelse(prob_test_lasso > 0.5, 1, 0)
                     
                     t2 = Sys.time()
                     time_lasso = as.numeric(difftime(t2,t1,units = "secs"))
                     
                     acc_train_lasso <- mean(pred_train_lasso == y_train)
                     acc_test_lasso  <- mean(pred_test_lasso == y_test)
                     
                     auc_train_lasso <- suppressMessages(auc(y_train, prob_train_lasso))
                     auc_test_lasso  <- suppressMessages(auc(y_test, prob_test_lasso))
                     
                     # ------------ swaglm
                     t3 = Sys.time()
                     alpha=0.15
                     pmax=10
                     fit_swaglm = swaglm(X_train, y_train, p_max = pmax, family = binomial(),
                                         alpha =alpha, verbose = FALSE)
                     y_pred_train_swaglm = predict(fit_swaglm, X_train)  
                     y_pred_test_swaglm  = predict(fit_swaglm, X_test)
                     
                     y_pred_train_swaglm_majority_class = ifelse(rowMeans(y_pred_train_swaglm$mat_reponse_prediction)>.5, 1, 0)
                     y_pred_test_swaglm_majority_class  = ifelse(rowMeans(y_pred_test_swaglm$mat_reponse_prediction)>.5, 1, 0)
                     
                     y_pred_train_swaglm_mean_pred = rowMeans(y_pred_train_swaglm$mat_reponse_prediction)
                     y_pred_test_swaglm_mean_pred  = rowMeans(y_pred_test_swaglm$mat_reponse_prediction)
                     
                     t4 = Sys.time()
                     time_swaglm = as.numeric(difftime(t4,t3,units = "secs"))
                     
                     acc_train_swaglm <- mean(y_pred_train_swaglm_majority_class == y_train)
                     acc_test_swaglm  <- mean(y_pred_test_swaglm_majority_class == y_test)
                     
                     auc_train_swaglm <- suppressMessages(auc(y_train, y_pred_train_swaglm_mean_pred))
                     auc_test_swaglm  <- suppressMessages(auc(y_test, y_pred_test_swaglm_mean_pred))
                     
                     # Return one row
                     c(seed = seed,
                       lambda = lambda_opt,
                       train_accuracy_lasso = acc_train_lasso,
                       test_accuracy_lasso  = acc_test_lasso,
                       train_auc_lasso      = auc_train_lasso,
                       test_auc_lasso       = auc_test_lasso,
                       time_lasso           = time_lasso,
                       train_accuracy_swaglm= acc_train_swaglm,
                       test_accuracy_swaglm = acc_test_swaglm,
                       train_auc_swaglm     = auc_train_swaglm,
                       test_auc_swaglm      = auc_test_swaglm,
                       time_swaglm          = time_swaglm)
                   }

# Stop cluster
stopCluster(cl)

# Store results as data.frame
mat_res <- as.data.frame(mat_res)
str(mat_res)



# make sound when finished
beepr::beep()

# plot boxplot of accuracy
boxplot(mat_res$train_accuracy_lasso, mat_res$train_accuracy_swaglm,
        mat_res$test_accuracy_lasso, mat_res$test_accuracy_swaglm,
        names= c("Lasso (train)", "SWAG (train)", "Lasso (test)", "SWAG (test)"), las = 1,
        ylab="Accuracy"
)
        

# plot boxplot of AUC 
boxplot(
        mat_res$train_auc_lasso, mat_res$train_auc_swaglm,
        mat_res$test_auc_lasso, mat_res$test_auc_swaglm,
        names= c("Lasso (train)", "SWAG (train)", "Lasso (test)", "SWAG (test)"), las = 1,
        ylab = "AUC"
)
     

# plot boxplot of time
boxplot(mat_res$time_lasso, mat_res$time_swaglm,
        names= c("Lasso", "SWAG"), las = 1
)


