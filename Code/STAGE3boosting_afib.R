## Jericho Lawson
## Summer 2019
## AFib Classification using Boosting/XGBoost (Used for parameter tuning for boosting).

# Pick specific version and covariates that you want to potentially include.
VERSION = "MIT"
#VERSION = "single_2017"

# Covariates will only be factored into PCA models.
COVARIATES = "nine"
#COVARIATES = "twelve"

# Pick the amount of folds or cross-validation.
k = 5; 
#k = 100
#k = dim(samp_data)[1]

#### Libraries ####
library(gbm); library(xgboost)

#### Data Collection ####

# Retrieves transition data from either the MIT-BIH dataset or the 2017 PhysioNet dataset.
if (VERSION == "MIT"){
  samp_data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_trans_dataMIT.csv", header = T)
}else if (VERSION == "single_2017"){
  samp_data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/Training_Set/all_tran_mat_2017single_2017.csv", header = T)
}

# Set-up for cross validation.
samp_data$State = as.factor(samp_data$State)
samp_data$State = factor(samp_data$State, levels = c('N', 'AFIB'))

set.seed(3); rand = sample(rep(1:k, length.out = dim(samp_data)[1]))

#### Functions ####

# Boosting testing.
boo_atr_funct = function(formula, train, test, ntrees, shrink, in_depth){
  start = Sys.time()
  train$State = ifelse(train$State == "AFIB", 1, 0)
  test$State = ifelse(test$State == "AFIB", 1, 0)
  
  set.seed(3)
  boo = gbm(formula, data = train, distribution = "bernoulli", n.trees = ntrees, shrinkage = shrink, interaction.depth = in_depth)
  boo_pred = predict(boo, test, n.trees=ntrees, type = "response")
  boo_pred = ifelse(boo_pred >= 0.5, "AFIB", "N")
  
  train$State = ifelse(train$State == 1, "AFIB", "N")
  test$State = ifelse(test$State == 1, "AFIB", "N")
  return(list(mean(boo_pred == test$State), test$State, boo_pred, Sys.time() - start))
}

# Function to run all boosting models.
atr_funct = function(i, rand, formula){
  samp = which(rand == i); train = samp_data[-samp, ]; test = samp_data[samp, ]
  
  #boo_res1 = boo_atr_funct(formula, train, test, 100, 0.001, 1)
  #boo_res2 = boo_atr_funct(formula, train, test, 5000, 0.001, 1)
  #boo_res3 = boo_atr_funct(formula, train, test, 100, 0.01, 1)
  #boo_res4 = boo_atr_funct(formula, train, test, 5000, 0.01, 1)
  #boo_res5 = boo_atr_funct(formula, train, test, 100, 0.01, 4)
  boo_res6 = boo_atr_funct(formula, train, test, 5000, 0.01, 4)
  #boo_res7 = boo_atr_funct(formula, train, test, 100, 0.001, 4)
  #boo_res8 = boo_atr_funct(formula, train, test, 5000, 0.001, 4)
  
  return(list(boo_res6))
}

# Function to run PCA boosting models.
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
  }
  
  cor_merged = prcomp(train_mat, retx = TRUE, center = TRUE, scale = TRUE); raw_merged = prcomp(train_mat, retx = TRUE, scale = FALSE) 
  train_pca = (train_mat) %*% raw_merged$rotation; test_pca = test_mat %*% raw_merged$rotation
  #cor_mat1 = cor(train_mat); lambdas = eigen(cor_mat1)$values; var_exp = cumsum(lambdas) / sum(lambdas)
  PCA = 5
  train = as.data.frame(train_pca[, 1:PCA]); test = as.data.frame(test_pca[, 1:PCA])
  train$State = samp_data[-samp, 5]; test$State = samp_data[samp, 5]
  
  #boo_res1 = boo_atr_funct(State ~ ., train, test, 100, 0.001, 1)
  #boo_res2 = boo_atr_funct(State ~ ., train, test, 5000, 0.001, 1)
  #boo_res3 = boo_atr_funct(State ~ ., train, test, 100, 0.01, 1)
  #boo_res4 = boo_atr_funct(State ~ ., train, test, 5000, 0.01, 1)
  #boo_res5 = boo_atr_funct(State ~ ., train, test, 100, 0.01, 4)
  boo_res6 = boo_atr_funct(State ~ ., train, test, 5000, 0.01, 4)
  #boo_res7 = boo_atr_funct(State ~ ., train, test, 100, 0.001, 4)
  #boo_res8 = boo_atr_funct(State ~ ., train, test, 5000, 0.001, 4)
  
  return(list(boo_res6))
}

#### Runs ####
boo_nine = lapply(1:5, atr_funct, rand = rand, formula = State ~ S.S + S.Reg + S.L + Reg.S + Reg.Reg + Reg.L + L.S + L.Reg + L.L + 
                    R.R.Length.Variance + R.R.Difference.Variance + R.R.Mean.Difference)
boo_rtl = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.L + R.R.Length.Variance)
boo_rtr = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.L + R.R.Length.Variance + S.Reg + R.R.Difference.Variance)
boo_pca5 = lapply(1:5, pca_funct, rand = rand)

#### Data organization. ####
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
    accs_nine = c(accs_nine, boo_nine[[fold]][[change]][[1]])
    preds_nine = c(preds_nine, boo_nine[[fold]][[change]][[2]])
    truths_nine = c(truths_nine, boo_nine[[fold]][[change]][[3]])
    times_nine = c(times_nine, boo_nine[[fold]][[change]][[4]])
    
    accs_rtl = c(accs_rtl, boo_rtl[[fold]][[change]][[1]])
    preds_rtl = c(preds_rtl, boo_rtl[[fold]][[change]][[2]])
    truths_rtl = c(truths_rtl, boo_rtl[[fold]][[change]][[3]])
    times_rtl = c(times_rtl, boo_rtl[[fold]][[change]][[4]])
    
    accs_rtr = c(accs_rtr, boo_rtr[[fold]][[change]][[1]])
    preds_rtr = c(preds_rtr, boo_rtr[[fold]][[change]][[2]])
    truths_rtr = c(truths_rtr, boo_rtr[[fold]][[change]][[3]])
    times_rtr = c(times_rtr, boo_rtr[[fold]][[change]][[4]])
    
    accs_pca5 = c(accs_pca5, boo_pca5[[fold]][[change]][[1]])
    preds_pca5 = c(preds_pca5, boo_pca5[[fold]][[change]][[2]])
    truths_pca5 = c(truths_pca5, boo_pca5[[fold]][[change]][[3]])
    times_pca5 = c(times_pca5, boo_pca5[[fold]][[change]][[4]])
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

#### Results ####

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