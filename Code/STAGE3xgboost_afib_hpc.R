## Jericho Lawson
## Summer 2019
## AFib Classification using XGBoost

# Pick specific version and covariates that you want to potentially include.
VERSION = "MIT"
#VERSION = "single_2017"

# Covariates will only be factored into PCA models.
COVARIATES = "nine"       # nine transition states
#COVARIATES = "twelve"      # nine transition states + R-R length variance + dRR mean + dRR variance

# Pick the amount of folds or cross-validation.
k = 5
#k = 100
#k = dim(samp_data)[1]

#### Libraries ####
library(gbm); library(xgboost)

#### Data Collection ####

# Retrieves transition data from either the MIT-BIH dataset or the 2017 PhysioNet dataset.
if (VERSION == "MIT"){
  samp_data = read.csv("/rhome/jlaws011/AFib/all_trans_dataMIT.csv", header = T)
}else if (VERSION == "single_2017"){
  samp_data = read.csv("/rhome/jlaws011/AFib/all_tran_mat_2017single_2017.csv", header = T)
}

samp_data$State = as.factor(samp_data$State)
samp_data$State = factor(samp_data$State, levels = c('N', 'AFIB'))
afib = samp_data[samp_data$State == "AFIB",]; nfib = samp_data[samp_data$State == "N",]

set.seed(3); rand = sample(rep(1:k, length.out = dim(samp_data)[1]))

#### Functions ####

# Runs XGBoost and classifies training data as AFib or non-AFib.
xgb_atr_funct = function(index, train, test, nrounds, eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, lambda, alpha){
  start = Sys.time()
  train$State = ifelse(train$State == "AFIB", 1, 0); train = data.matrix(train)
  test$State = ifelse(test$State == "AFIB", 1, 0); test = data.matrix(test)
  
  set.seed(3)
  xgb = xgboost(data = as.matrix(train[, index]), label = train[,5], nrounds = nrounds, eta = eta, gamma = gamma,
                max_depth = max_depth, min_child_weight = min_child_weight, subsample = subsample,
                colsample_bytree = colsample_bytree, lambda = lambda, alpha = alpha, 
                scale_pos_weight = sum(train[, 5] == 0) / sum(train[, 5] == 1),
                objective = "binary:logistic", verbose = 0)
  xgb_pred = predict(xgb, as.matrix(test[, index]))
  xgb_pred = ifelse(xgb_pred >= 0.5, "AFIB", "N")
  
  test[,5] = ifelse(test[,5] == 1, "AFIB", "N")
  
  return(list(mean(xgb_pred == test[,5]), test[,5], xgb_pred, Sys.time() - start))
}

# Runs XGBoost and classifies PCA dimensions from training data as AFib or non-AFib.
xgb_pca_funct = function(train, test, nrounds, eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, lambda, alpha, dims){
  start = Sys.time()
  train$State = ifelse(train$State == "AFIB", 1, 0); train = data.matrix(train)
  test$State = ifelse(test$State == "AFIB", 1, 0); test = data.matrix(test)
  
  set.seed(3)
  xgb = xgboost(data = train[, 1:dims], label = train[,(dims + 1)], nrounds = nrounds, eta = eta, gamma = gamma,
                max_depth = max_depth, min_child_weight = min_child_weight, subsample = subsample,
                colsample_bytree = colsample_bytree, lambda = lambda, alpha = alpha, 
                scale_pos_weight = sum(train[, (dims + 1)] == 0) / sum(train[, (dims + 1)] == 1),
                objective = "binary:logistic", verbose = 0)
  xgb_pred = predict(xgb, test[, 1:dims])
  xgb_pred = ifelse(xgb_pred >= 0.5, "AFIB", "N")
  
  test[,(dims + 1)] = ifelse(test[,(dims + 1)] == 1, "AFIB", "N")
  
  return(list(mean(xgb_pred == test[,(dims + 1)]), test[,(dims + 1)], xgb_pred, Sys.time() - start))
}

# Function to run all XGBoost models.
atr_funct = function(i, rand, index){
  
  samp = which(rand == i); train = samp_data[-samp, ]; test = samp_data[samp, ]
  
  xgb_res1 = xgb_atr_funct(index, train, test, 5000, 0.01, 0, 4, 1, 1, 1, 0, 1)
  #xgb_res2 = xgb_atr_funct(index, train, test, 100, 0.3, 0, 4, 1, 1, 1, 0, 1)
  #xgb_res3 = xgb_atr_funct(index, train, test, 100, 0.3, 0, 6, 1, 1, 1, 0, 1)
  #xgb_res4 = xgb_atr_funct(index, train, test, 100, 0.3, 0, 8, 1, 1, 1, 0, 1)
  #xgb_res5 = xgb_atr_funct(index, train, test, 250, 0.3, 0, 2, 1, 1, 1, 0, 1)
  #xgb_res6 = xgb_atr_funct(index, train, test, 250, 0.3, 0, 4, 1, 1, 1, 0, 1)
  #xgb_res7 = xgb_atr_funct(index, train, test, 250, 0.3, 0, 6, 1, 1, 1, 0, 1)
  #xgb_res8 = xgb_atr_funct(index, train, test, 250, 0.3, 0, 8, 1, 1, 1, 0, 1)
  
  #return(list(xgb_res1, xgb_res2, xgb_res3, xgb_res4, xgb_res5, xgb_res6, xgb_res7, xgb_res8))
  return(list(xgb_res1))
}

# Function run all XGBoost models using PCA.
pca_funct = function(i, rand){
  samp = which(rand == i)
  
  if (COVARIATES == "nine"){
    train_mat = as.matrix(samp_data[-samp, 8:16]); test_mat = as.matrix(samp_data[samp, 8:16])
  }else if (COVARIATES == "twelve"){
    if (VERSION == "MIT"){
      train_mat = as.matrix(samp_data[-samp, 8:19]); test_mat = as.matrix(samp_data[samp, 8:19])
    }else{
      train_mat = as.matrix(samp_data[-samp, c(8:16, 18:20)]); test_mat = as.matrix(samp_data[samp, c(8:16, 18:20)])
    }
  }else{
    train_mat = as.matrix(samp_data[-samp, 8:22]); test_mat = as.matrix(samp_data[samp, 8:22]) 
  }
  
  cor_merged = prcomp(train_mat, retx = TRUE, center = TRUE, scale = TRUE); raw_merged = prcomp(train_mat, retx = TRUE, scale = FALSE) 
  train_pca = (train_mat) %*% raw_merged$rotation; test_pca = test_mat %*% raw_merged$rotation
  #cor_mat1 = cor(train_mat); lambdas = eigen(cor_mat1)$values; var_exp = cumsum(lambdas) / sum(lambdas)
  PCA = 5
  train = as.data.frame(train_pca[, 1:PCA]); test = as.data.frame(test_pca[, 1:PCA])
  train$State = samp_data[-samp, 5]; test$State = samp_data[samp, 5]
  
  xgb_res1 = xgb_pca_funct(train, test, 5000, 0.01, 0, 4, 1, 1, 1, 0, 1, PCA)
  #xgb_res2 = xgb_pca_funct(train, test, 100, 0.3, 0, 4, 1, 1, 1, 0, 1, PCA)
  #xgb_res3 = xgb_pca_funct(train, test, 100, 0.3, 0, 6, 1, 1, 1, 0, 1, PCA)
  #xgb_res4 = xgb_pca_funct(train, test, 100, 0.3, 0, 8, 1, 1, 1, 0, 1, PCA)
  #xgb_res5 = xgb_pca_funct(train, test, 250, 0.3, 0, 2, 1, 1, 1, 0, 1, PCA)
  #xgb_res6 = xgb_pca_funct(train, test, 250, 0.3, 0, 4, 1, 1, 1, 0, 1, PCA)
  #xgb_res7 = xgb_pca_funct(train, test, 250, 0.3, 0, 6, 1, 1, 1, 0, 1, PCA)
  #xgb_res8 = xgb_pca_funct(train, test, 250, 0.3, 0, 8, 1, 1, 1, 0, 1, PCA)
  
  #return(list(xgb_res1, xgb_res2, xgb_res3, xgb_res4, xgb_res5, xgb_res6, xgb_res7, xgb_res8))
  return(list(xgb_res1))
}

#### Runs ####
if (COVARIATES == "nine"){
  xgb_nine = lapply(1:k, atr_funct, rand = rand, index = 8:16) 
  xgb_rtl = lapply(1:k, atr_funct, rand = rand, index = c(13))
  xgb_rtr = lapply(1:k, atr_funct, rand = rand, index = c(9, 13))
  xgb_pca5 = lapply(1:k, pca_funct, rand = rand)
}else if (COVARIATES == "twelve"){
  if (VERSION == "MIT"){
    xgb_nine = lapply(1:k, atr_funct, rand = rand, index = 8:19) 
    xgb_rtl = lapply(1:k, atr_funct, rand = rand, index = c(13, 17))
    xgb_rtr = lapply(1:k, atr_funct, rand = rand, index = c(9, 13, 17, 18))
    xgb_pca5 = lapply(1:k, pca_funct, rand = rand)
  }else{
    xgb_nine = lapply(1:k, atr_funct, rand = rand, index = c(8:16, 18:20)) 
    xgb_rtl = lapply(1:k, atr_funct, rand = rand, index = c(13, 18))
    xgb_rtr = lapply(1:k, atr_funct, rand = rand, index = c(9, 13, 18, 19))
    xgb_pca5 = lapply(1:k, pca_funct, rand = rand)
  }
}

#### Data organization ####
all_accs_nine = c(); all_sens_nine = c(); all_specs_nine = c(); all_sds_nine = c(); all_f1_nine = c(); all_times_nine = c()
all_accs_rtl = c(); all_sens_rtl = c(); all_specs_rtl = c(); all_sds_rtl = c(); all_f1_rtl = c(); all_times_rtl = c()
all_accs_rtr = c(); all_sens_rtr = c(); all_specs_rtr = c(); all_sds_rtr = c(); all_f1_rtr = c(); all_times_rtr = c()
all_accs_pca5 = c(); all_sens_pca5 = c(); all_specs_pca5 = c(); all_sds_pca5 = c(); all_f1_pca5 = c(); all_times_pca5 = c()

for (change in 1:1){
  accs_nine = c(); preds_nine = c(); truths_nine = c(); times_nine = c()
  accs_rtl = c(); preds_rtl = c(); truths_rtl = c(); times_rtl = c()
  accs_rtr = c(); preds_rtr = c(); truths_rtr = c(); times_rtr = c()
  accs_pca5 = c(); preds_pca5 = c(); truths_pca5 = c(); times_pca5 = c()
  for (fold in 1:5){
    accs_nine = c(accs_nine, xgb_nine[[fold]][[change]][[1]])
    preds_nine = c(preds_nine, xgb_nine[[fold]][[change]][[3]])
    truths_nine = c(truths_nine, xgb_nine[[fold]][[change]][[2]])
    times_nine = c(times_nine, xgb_nine[[fold]][[change]][[4]])
    
    accs_rtl = c(accs_rtl, xgb_rtl[[fold]][[change]][[1]])
    preds_rtl = c(preds_rtl, xgb_rtl[[fold]][[change]][[3]])
    truths_rtl = c(truths_rtl, xgb_rtl[[fold]][[change]][[2]])
    times_rtl = c(times_rtl, xgb_rtl[[fold]][[change]][[4]])
    
    accs_rtr = c(accs_rtr, xgb_rtr[[fold]][[change]][[1]])
    preds_rtr = c(preds_rtr, xgb_rtr[[fold]][[change]][[3]])
    truths_rtr = c(truths_rtr, xgb_rtr[[fold]][[change]][[2]])
    times_rtr = c(times_rtr, xgb_rtr[[fold]][[change]][[4]])
    
    accs_pca5 = c(accs_pca5, xgb_pca5[[fold]][[change]][[1]])
    preds_pca5 = c(preds_pca5, xgb_pca5[[fold]][[change]][[3]])
    truths_pca5 = c(truths_pca5, xgb_pca5[[fold]][[change]][[2]])
    times_pca5 = c(times_pca5, xgb_pca5[[fold]][[change]][[4]])
  }
  all_accs_nine = c(all_accs_nine, mean(accs_nine))
  all_sds_nine = c(all_sds_nine, sd(accs_nine))
  tab_nine = table(preds_nine, truths_nine)
  all_sens_nine = c(all_sens_nine, tab_nine[1, 1] / sum(tab_nine[, 1]))
  all_specs_nine = c(all_specs_nine, tab_nine[2, 2] / sum(tab_nine[, 2]))
  all_f1_nine = c(all_f1_nine, 2 * ((tab_nine[1, 1] / sum(tab_nine[1,])) * (tab_nine[1, 1] / sum(tab_nine[,1]))) / 
                    (((tab_nine[1, 1] / sum(tab_nine[1,])) + (tab_nine[1, 1] / sum(tab_nine[,1])))))
  all_times_nine = c(all_times_nine, sum(times_nine))
  
  all_accs_rtl = c(all_accs_rtl, mean(accs_rtl))
  all_sds_rtl = c(all_sds_rtl, sd(accs_rtl))
  tab_rtl = table(preds_rtl, truths_rtl)
  all_sens_rtl = c(all_sens_rtl, tab_rtl[1, 1] / sum(tab_rtl[, 1]))
  all_specs_rtl = c(all_specs_rtl, tab_rtl[2, 2] / sum(tab_rtl[, 2]))
  all_f1_rtl = c(all_f1_rtl, 2 * ((tab_rtl[1, 1] / sum(tab_rtl[1,])) * (tab_rtl[1, 1] / sum(tab_rtl[,1]))) / 
                   (((tab_rtl[1, 1] / sum(tab_rtl[1,])) + (tab_rtl[1, 1] / sum(tab_rtl[,1])))))
  all_times_rtl = c(all_times_rtl, sum(times_rtl))
  
  all_accs_rtr = c(all_accs_rtr, mean(accs_rtr))
  all_sds_rtr = c(all_sds_rtr, sd(accs_rtr))
  tab_rtr = table(preds_rtr, truths_rtr)
  all_sens_rtr = c(all_sens_rtr, tab_rtr[1, 1] / sum(tab_rtr[, 1]))
  all_specs_rtr = c(all_specs_rtr, tab_rtr[2, 2] / sum(tab_rtr[, 2]))
  all_f1_rtr = c(all_f1_rtr, 2 * ((tab_rtr[1, 1] / sum(tab_rtr[1,])) * (tab_rtr[1, 1] / sum(tab_rtr[,1]))) / 
                   (((tab_rtr[1, 1] / sum(tab_rtr[1,])) + (tab_rtr[1, 1] / sum(tab_rtr[,1])))))
  all_times_rtr = c(all_times_rtr, sum(times_rtr))
  
  all_accs_pca5 = c(all_accs_pca5, mean(accs_pca5))
  all_sds_pca5 = c(all_sds_pca5, sd(accs_pca5))
  tab_pca5 = table(preds_pca5, truths_pca5)
  all_sens_pca5 = c(all_sens_pca5, tab_pca5[1, 1] / sum(tab_pca5[, 1]))
  all_specs_pca5 = c(all_specs_pca5, tab_pca5[2, 2] / sum(tab_pca5[, 2]))
  all_f1_pca5 = c(all_f1_pca5, 2 * ((tab_pca5[1, 1] / sum(tab_pca5[1,])) * (tab_pca5[1, 1] / sum(tab_pca5[,1]))) / 
                    (((tab_pca5[1, 1] / sum(tab_pca5[1,])) + (tab_pca5[1, 1] / sum(tab_pca5[,1])))))
  all_times_pca5 = c(all_times_pca5, sum(times_pca5))
}

#### Results (printed out) ####
paste("Model #1: Accuracies:", all_accs_nine)
paste("Model #1: SE of Accuracies:", all_sds_nine)
paste("Model #1: Sensitivities:", all_sens_nine)
paste("Model #1: Specificities:", all_specs_nine)
paste("Model #1: F-1 Scores:", all_f1_nine)
paste("Model #1: Comp. Times:", all_times_nine)

paste("Model #2: Accuracies:", all_accs_rtl)
paste("Model #2: SE of Accuracies:", all_sds_rtl)
paste("Model #2: Sensitivities:", all_sens_rtl)
paste("Model #2: Specificities:", all_specs_rtl)
paste("Model #2: F-1 Scores:", all_f1_rtl)
paste("Model #2: Comp. Times:", all_times_rtl)

paste("Model #3: Accuracies:", all_accs_rtr)
paste("Model #3: SE of Accuracies:", all_sds_rtr)
paste("Model #3: Sensitivities:", all_sens_rtr)
paste("Model #3: Specificities:", all_specs_rtr)
paste("Model #3: F-1 Scores:", all_f1_rtr)
paste("Model #3: Comp. Times:", all_times_rtr)

paste("Model #4: Accuracies:", all_accs_pca5)
paste("Model #4: SE of Accuracies:", all_sds_pca5)
paste("Model #4: Sensitivities:", all_sens_pca5)
paste("Model #4: Specificities:", all_specs_pca5)
paste("Model #4: F-1 Scores:", all_f1_pca5)
paste("Model #4: Comp. Times:", all_times_pca5)

###################################################################################################
