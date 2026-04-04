library("VariantAnnotation")
library("tidyverse")
library("ggplot2")
library("ggfortify")
library(glmnet)
library(caret)
library(pROC)
library(readxl)
library(dplyr)
library("stringr")
library(rpart)
library(randomForest)
library(rpart.plot)
library(rattle)
library(ipred)
#Load data in the file
SLC24A5<- readVcf("gene vcf files/SLC24A5.vcf.gz")
EDAR <- readVcf("gene vcf files/EDAR.vcf.gz")
#Data Exploration
SLC24A5
samples(header(SLC24A5))
seqlevels(rowRanges(SLC24A5))
rowRanges(SLC24A5)
#Connvert to data frame
DF_SLC24A5 <-  as.data.frame(geno(SLC24A5)$GT)
DF_EDAR <- as.data.frame(geno(EDAR)$GT)

#Function to convert genotype information to numeric allele dosage
#Homozygous major allele 0, Heterozygous 1, Homozygous minor allele 2.
convert_gt <- function(x) {
  x <- as.character(x)
  x[x %in% c("0|0", "0/0")] <- "0"
  x[x %in% c("0|1", "1|0", "0/1", "1/0")] <- "1"
  x[x %in% c("1|1", "1/1")] <- "2"
  x[x %in% c("./.", ".|.")] <- NA
  as.numeric(x)
}
#Applying function 
DF_SLC24A5[] <- as.data.frame(lapply(DF_SLC24A5, convert_gt))
DF_EDAR[] <- as.data.frame(lapply(DF_EDAR, convert_gt))

#Adding metadata 
metadata <- read.table(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",
  header = TRUE,
  stringsAsFactors = FALSE
)

head(metadata)

#Transpose SNP data to have sample IDs as rownames
DF_SLC24A5 <- as.data.frame(t(DF_SLC24A5))
DF_EDAR <- as.data.frame(t(DF_EDAR))


#Combine all data frames
Genotype_df <-  cbind(metadata,DF_SLC24A5,DF_EDAR)

#Working with the European  an East asian superpopulation (AFR)
#We intend to use SNP data from these Ancestry informative genes to classify the population.
COMBINED_DF <- Genotype_df %>%
  filter(super_pop %in% c("EAS", "EUR"))
#Data exploration
ggplot(COMBINED_DF, aes(x = super_pop)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Populations within dataset", x = "Superpopulation", y = "Count")

#Remove SNPs with NA values
Cleaned_DF <- COMBINED_DF[, colSums(is.na(COMBINED_DF)) == 0, drop = FALSE]


# Remove columns(SNPs) with zero variance
X <- Cleaned_DF[, -(1:4)]

X_filtered <- X[, apply(X, 2, var) != 0]

#Filtering low minor allele frequency below 0.01
maf <- apply(X_filtered, 2, function(snp) {
  p <- sum(snp, na.rm = TRUE) / (2 * sum(!is.na(snp)))
  min(p, 1 - p)
})

X_filtered <- X_filtered[, maf > 0.01]


# Run PCA
PCA <- prcomp(X_filtered, center = TRUE)
var_explained <- (PCA$sdev^2) / sum(PCA$sdev^2)
percent_var <- var_explained * 100

percent_var
#Combine data with metadata for plot
metadata_sub <-  metadata %>%
  filter(super_pop %in% c("EAS", "EUR"))

X_filtered2 <- cbind(metadata_sub, X_filtered)

#Plot
autoplot( PCA , data = X_filtered2 , colour = "super_pop", main = "PCA: PC1 vs PC2" )

#Method 1
#Logistic regression and cross validation
## Read and prepare data----


Xd <- model.matrix(~ . - 1, data = X_filtered)
Y <- as.factor(X_filtered2$super_pop)



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


# Cross-Validation----
# Define alpha search grid
alpha_grid <- seq(0, 1, .01)

# Create data frames for collecting cross-validation results
result_10_fold <- data.frame()
result_20_fold <- data.frame()


# 10-Fold
for (a in alpha_grid) {
  cv_10fold <- cv.glmnet(X_train, Y_train, 
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
  labs(title = "Deviance vs Alpha, 10-Fold",
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

# Extract top 11 predictors for each model. 
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


#Method 2 
#Random forest with mtry
#Using X_filtered2 our Df with metadata
set.seed(4242) 
#Create 70/30 training test split using an index
train.index <- sample(1:nrow(X_filtered2), round(0.70*nrow(X_filtered2),0))

cols_to_drop <- c("sample","pop","gender")

X_clean_RF <- X_filtered2[, !(names(X_filtered2) %in% cols_to_drop)]

train_data <- X_clean_RF[train.index, ]
test_data  <- X_clean_RF[-train.index, ]

#Tune random forest using mtry
#mtry set to 10,20,30,50 due to number of predictors

set.seed(2445)


control <- trainControl(method="repeatedcv", 
                        number=5, 
                        repeats=3,
                        allowParallel = FALSE)


modelRF <- train(
  super_pop ~ . ,
  data = train_data,
  method = "rf",
  trControl = control,
  tuneGrid = expand.grid(mtry = c(10, 20, 30, 50))
)
plot(modelRF)

#Variable importance
imp <- as.data.frame(varImp(modelRF)$importance)

imp$SNP <- rownames(imp)
imp$Overall <- rowMeans(imp[, sapply(imp, is.numeric)])
#Top 20 SNPs
top20 <- imp[order(imp$Overall, decreasing = TRUE), ]
top20 <- head(top20, 20)
top20
#plot
ggplot(top20, aes(x = reorder(SNP, Overall), y = Overall)) +
  geom_col(fill = 'brown') +
  coord_flip() +
  labs(x = "SNP", y = "Variable Importance", title = "Top 20 SNPs from random forest") + 
  theme_minimal()

#Predict
RFpredict <- predict(
  modelRF,
  newdata = test_data,
  type = "prob"
)[,1]

roc_obj <- roc(X_filtered2[-train.index,]$super_pop,RFpredict)
roc_obj
auc(roc_obj)
#Perfect classifier due to highly informative SNPs.

#Error rate 
#Training error
tail(modelRF$finalModel$err.rate, 1)

#Test error
pred_class <- predict(modelRF, newdata = test_data)

mean(pred_class != test_data$super_pop)

