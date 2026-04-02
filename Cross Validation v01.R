library(glmnet)
library(ggplot2)
library(caret)
library(pROC)
library(readxl)
library(dplyr)
library("stringr")
library(tidyverse)
library(ggfortify)


#Problem 2
## Read and prepare data----


Xd <- model.matrix(~ . - 1, data = X_filt_head)
Y <- as.factor(Y_head$pop)



# Set seed for cross-validation (for sanity purposes)
set.seed(1717)

# Split data into training and test data (75/25)

train_index <- createDataPartition(Y, p = 0.75, list = FALSE)
X_train_prescaled <- Xd[train_index, ]
X_test_prescaled <- Xd[-train_index, ]
Y_train <- Y[train_index]
Y_test <- Y[-train_index]


# Standardize predictors (glmnet does this by default with standardize=TRUE)
# But we'll scale manually for consistency across methods
X_train <- scale(X_train_prescaled)
X_test <- scale(X_test_prescaled, 
                center = attr(X_train, "scaled:center"),
                scale = attr(X_train, "scaled:scale"))

# Remove any columns with NA values (if scaling creates them due to zero variance)
na_cols <- which(colSums(is.na(X_train)) > 0)
if (length(na_cols) > 0) {
  X_train_scaled <- X_train[, -na_cols]
  X_test_scaled <- X_test[, -na_cols]
  cat("Removed", length(na_cols), "columns with zero variance\n")
}

# Cross-Validation----

# NOTE: I used "deviance" instead of "auc" because the R warning message said too few observations per fold when using "auc".  Try out both methods, though, if you wish.

# Define alpha search grid
alpha_grid <- seq(0, 1, .01)

# Create data frames for collecting cross-validation results
result_10_fold <- data.frame()
result_20_fold <- data.frame()


# 10-Fold
for (a in alpha_grid) {
  cv_10fold <- cv.glmnet(X_train_scaled, Y_train, 
                         alpha = a, 
                         nfolds = 10, 
                         family = "multinomial", 
                         type.measure = "deviance", 
                         keep = TRUE)
  
  result_10_fold <- rbind(result_10_fold, data.frame(
    alpha = a, 
    lambda_min = cv_10fold$lambda.min, 
    lambda_1se = cv_10fold$lambda.1se, 
    cvm_min = min(cv_10fold$cvm), 
    cvm_1se = cv_10fold$cvm[which(cv_10fold$lambda == cv_10fold$lambda.1se)],
    #adding in nzero, as was included in original for loop
    nzero_min = cv_10fold$nzero[which(cv_10fold$lambda == cv_10fold$lambda.min)],
    nzero_1se = cv_10fold$nzero[which(cv_10fold$lambda == cv_10fold$lambda.1se)]
  ))
}

print(result_10_fold)

# Extract best alpha and lambda values (10-Fold) by min CV error
best_alpha_min_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_min)]
best_alpha_1se_10 <- result_10_fold$alpha[which.min(result_10_fold$cvm_1se)]

best_lambda_min_10 <- result_10_fold$lambda_min[result_10_fold$alpha == best_alpha_min_10]
best_lambda_1se_10 <- result_10_fold$lambda_1se[result_10_fold$alpha == best_alpha_1se_10]


#Plots for visualization and alpha grid search
plot(cv_10fold, main = "Elastic Net (10-Fold) - Deviance")


ggplot(result_10_fold, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "CV Error vs Alpha, 10-Fold",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "Deviance",
       color = "Lambda Type") +
  theme_minimal()

# 20-Fold
for (a in alpha_grid) {
  cv_20fold <- cv.glmnet(X_train, Y_train, 
                         alpha = a, 
                         nfolds = 20, 
                         family = "binomial", 
                         type.measure = "deviance", 
                         keep = TRUE)
  
  result_20_fold <- rbind(result_20_fold, data.frame(
    alpha = a, 
    lambda_min = cv_20fold$lambda.min, 
    lambda_1se = cv_20fold$lambda.1se, 
    cvm_min = min(cv_20fold$cvm), 
    cvm_1se = cv_20fold$cvm[which(cv_20fold$lambda == cv_20fold$lambda.1se)],
    #adding in nzero, as was included in original for loop
    nzero_min = cv_20fold$nzero[which(cv_20fold$lambda == cv_20fold$lambda.min)],
    nzero_1se = cv_20fold$nzero[which(cv_20fold$lambda == cv_20fold$lambda.1se)]
  ))
}


print(result_20_fold)

# Extract best alpha and lambda values (20-Fold) by min CV error
best_alpha_min_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_min)]
best_alpha_1se_20 <- result_20_fold$alpha[which.min(result_20_fold$cvm_1se)]

best_lambda_min_20 <- result_20_fold$lambda_min[result_20_fold$alpha == best_alpha_min_20]
best_lambda_1se_20 <- result_20_fold$lambda_1se[result_20_fold$alpha == best_alpha_1se_20]

#Plots for visualization and alpha grid search
plot(cv_20fold, main = "Elastic Net (20-Fold) - Deviance")

ggplot(result_20_fold, aes(x = alpha)) +
  geom_line(aes(y = cvm_min, color = "lambda.min"), size = 1.2) +
  geom_line(aes(y = cvm_1se, color = "lambda.1se"), size = 1.2) +
  geom_point(aes(y = cvm_min, color = "lambda.min"), size = 3) +
  geom_point(aes(y = cvm_1se, color = "lambda.1se"), size = 3) +
  labs(title = "Deviance vs Alpha, 20-Fold",
       x = "Alpha (0=Ridge, 1=Lasso)",
       y = "Deviance",
       color = "Lambda Type") +
  theme_minimal()



# Fitting the final Elastic-Net models based on training data

# 10-Fold
model_10fold_min <- glmnet(X_train, Y_train, 
                           alpha = best_alpha_min_10, 
                           lambda = best_lambda_min_10,
                           lambda_type = "min",
                           family = "binomial")



model_10fold_1se <- glmnet(X_train, Y_train, 
                           alpha = best_alpha_1se_10, 
                           lambda = best_lambda_1se_10,
                           lambda_type = "1se",
                           family = "binomial")

# 20-Fold
model_20fold_min <- glmnet(X_train, Y_train, 
                           alpha = best_alpha_min_20, 
                           lambda = best_lambda_min_20,
                           lambda_type = "min",
                           family = "binomial")


model_20fold_1se <- glmnet(X_train, Y_train, 
                           alpha = best_alpha_1se_20, 
                           lambda = best_lambda_1se_20,
                           lambda_type = "1se",
                           family = "binomial")


# Function to extract top N predictors from a model
extract_top_predictors <- function(model, n = 10, lambda_type = "min") {
  # Get coefficients
  if (lambda_type == "min") {
    coefs <- as.matrix(coef(model))[, 1]
  } else {
    # For models stored differently
    coefs <- as.matrix(coef(model, s = paste0("lambda.", lambda_type)))[, 1]
  }
  
  # Remove intercept
  coefs_no_intercept <- coefs[-1]
  
  # Get absolute values for ranking
  abs_coefs <- abs(coefs_no_intercept)
  
  # Get top N predictors
  top_indices <- order(abs_coefs, decreasing = TRUE)[1:min(n, length(abs_coefs))]
  
  # Create result data frame
  result <- data.frame(
    Predictor = names(coefs_no_intercept)[top_indices],
    Coefficient = coefs_no_intercept[top_indices],
    Absolute_Effect = abs_coefs[top_indices],
    Rank = 1:length(top_indices)
  )
  
  return(result)
}

# Extract top 11 predictors for each model. 11 selected to try and identify impact of both age and the top 10 biomarkers (if that many)
m10_top_min <- extract_top_predictors(model_10fold_min, n = 11, lambda_type = "min")
print(m10_top_min)

m10_top_lse <- extract_top_predictors(model_10fold_1se, n = 11, lambda_type = "1se")
print(m10_top_lse)

m20_top_min <- extract_top_predictors(model_20fold_min, n = 11, lambda_type = "min")
print(m20_top_min)

m20_top_lse <- extract_top_predictors(model_20fold_1se, n = 11, lambda_type = "1se")
print(m20_top_lse)


# Function to calculate performance metrics
calculate_metrics <- function(model, X_test, Y_test, lambda_type = "min") {
  # Make predictions
  if (lambda_type == "min") {
    predictions <- predict(model, newx = X_test)
  } else {
    predictions <- predict(model, newx = X_test, s = paste0("lambda.", lambda_type))
  }
  
  # Calculate metrics
  mse <- mean((Y_test - predictions)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(Y_test - predictions))
  r_squared <- 1 - sum((Y_test - predictions)^2) / sum((Y_test - mean(Y_test))^2)
  
  # Count non-zero coefficients (excluding intercept)
  if (lambda_type == "min") {
    n_coefs <- sum(coef(model)[-1] != 0)
  } else {
    n_coefs <- sum(coef(model, s = paste0("lambda.", lambda_type))[-1] != 0)
  }
  
  return(c(MSE = mse, RMSE = rmse, MAE = mae, R2 = r_squared, NonZero = n_coefs))
}


# Calculate metrics for all models. We attempted to get a comparison table working for further cross-model metrics, but we could only get Non-Zero coefficient comparison working before the deadline.

comparison <- data.frame()

# 10-Fold
comparison <- rbind(comparison, 
                    data.frame(Model = "10-Fold (lambda.min)", 
                               as.list(calculate_metrics(model_10fold_min, X_test, Y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "10-Fold (lambda.1se)", 
                               as.list(calculate_metrics(model_10fold_1se, X_test, Y_test, "1se"))))

# 20-Fold
comparison <- rbind(comparison, 
                    data.frame(Model = "20-Fold (lambda.min)", 
                               as.list(calculate_metrics(model_20fold_min, X_test, Y_test, "min"))))

comparison <- rbind(comparison, 
                    data.frame(Model = "20-Fold (lambda.1se)", 
                               as.list(calculate_metrics(model_20fold_1se, X_test, Y_test, "1se"))))

print(comparison)
