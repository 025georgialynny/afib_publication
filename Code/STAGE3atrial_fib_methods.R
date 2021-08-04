## Jericho Lawson
## Summer 2019
## Atrial Fib - Resampling and Methods (Only used for MIT-BIH data and binary 2017 PhysioNet data.)

# Pick specific version and covariates that you want to potentially include.
VERSION = "MIT"
#VERSION = "single_2017"

COVARIATES = "nine"
#COVARIATES = "twelve"

# Pick the amount of folds or cross-validation.
k = 5; 
#k = 100
#k = dim(samp_data)[1]

#### Libraries ####
library(MASS); library(gbm); library(scales); library(ggplot2); library(reshape2); library(ggpubr); library(tidyr)

#### Data Collection ####

# Retrieves transition data from either the MIT-BIH dataset or the 2017 PhysioNet dataset.
if (VERSION == "MIT"){
  samp_data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_trans_dataMIT.csv", header = T)
}else if (VERSION == "single_2017"){
  samp_data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/Training_Set/all_tran_mat_2017single_2017.csv", header = T)
}

#### Set-up for Cross Validation ####

samp_data$State = as.factor(samp_data$State)
set.seed(3); rand = sample(rep(1:k, length.out = dim(samp_data)[1]))

#### Functions ####

# Logistic regression testing.
glm_atr_funct = function(formula, train, test){
  glm = glm(formula, data = train, family = "binomial")
  glm_probs = predict(glm, test, type = "response")
  glm_pred = rep("AFIB", dim(test)[1]); glm_pred[glm_probs > .5] = "N"
  return(list(mean(glm_pred == test$State), test$State, glm_pred))
}

# Linear discriminant analysis testing.
lda_atr_funct = function(formula, train, test){
  lda = lda(formula, data = train)
  lda_pred = predict(lda, test, type = "response")$class
  return(list(mean(lda_pred == test$State), test$State, lda_pred))
}

# Quadratic discriminant analysis testing.
qda_atr_funct = function(formula, train, test){
  qda = qda(formula, data = train)
  qda_pred = predict(qda, test, type = "response")$class
  return(list(mean(qda_pred == test$State), test$State, qda_pred))
}

# Boosting testing.
boo_atr_funct = function(formula, train, test){
  train$State = ifelse(train$State == "AFIB", 1, 0); test$State = ifelse(test$State == "AFIB", 1, 0)
  
  set.seed(3)
  boo = gbm(formula, data = train, distribution = "bernoulli", n.trees = 5000, shrinkage = 0.01, interaction.depth = 4)
  boo_pred = predict(boo, test, n.trees = 5000, type = "response")
  boo_pred = ifelse(boo_pred >= 0.5, "AFIB", "N")

  train$State = ifelse(train$State == 1, "AFIB", "N"); test$State = ifelse(test$State == 1, "AFIB", "N")
  return(list(mean(boo_pred == test$State), test$State, boo_pred))
}

# Function that runs all testing.
atr_funct = function(i, rand, formula){
  samp = which(rand == i); train = samp_data[-samp, ]; test = samp_data[samp, ]
  
  glm_res = glm_atr_funct(formula, train, test)
  lda_res = lda_atr_funct(formula, train, test)
  qda_res = qda_atr_funct(formula, train, test)
  boo_res = boo_atr_funct(formula, train, test)
  
  return(list(glm_res, lda_res, qda_res, boo_res))
}

# Produces accuracies, predictions, and true values from testing.
stats_funct = function(num, data_list){
  accs = c(); truths = c(); predicts = c()
  for (i in 1:length(data_list)){
    accs = c(accs, data_list[[i]][[num]][[1]])
    truths = c(truths, data_list[[i]][[num]][[2]])
    predicts = c(predicts, data_list[[i]][[num]][[3]])
  }
  con_tab = table(predicts, truths)
  return(list(mean(accs), sd(accs), con_tab, con_tab[1,1] / sum(con_tab[,1]), con_tab[2,2] / sum(con_tab[,2])))
}

#### Models and Results ####

if (COVARIATES == "nine"){
  nine_cov = lapply(1:5, atr_funct, rand = rand, formula = State ~ S.S + S.Reg + S.L + Reg.S + Reg.Reg + Reg.L + L.S + L.Reg + L.L)
  reg_long = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.L)
  reg_reg = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.Reg)
  
}else if (COVARIATES == "twelve"){
  nine_cov = lapply(1:5, atr_funct, rand = rand, formula = State ~ S.S + S.Reg + S.L + Reg.S + Reg.Reg + Reg.L + L.S + L.Reg + L.L + R.R.Length.Variance + R.R.Difference.Variance + R.R.Mean.Difference)
  reg_long = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.L + R.R.Length.Variance + R.R.Difference.Variance + R.R.Mean.Difference)
  reg_reg = lapply(1:5, atr_funct, rand = rand, formula = State ~ Reg.Reg + R.R.Length.Variance + R.R.Difference.Variance + R.R.Mean.Difference)
}


# Model 1: Nine Covariates
nine_cov_results = lapply(1:4, stats_funct, data_list = nine_cov)

# Model 2: Regular-to-Long Covariate
reg_long_results = lapply(1:4, stats_funct, data_list = reg_long)

# Model 3: Regular-to-Regular Covariate
reg_reg_results = lapply(1:4, stats_funct, data_list = reg_reg)

# Model 4: PCA Covariates
lda_diff_pca_accs = c(); lda_diff_pca_sd = c(); lda_diff_pca_sens = c(); lda_diff_pca_spec = c()
glm_diff_pca_accs = c(); glm_diff_pca_sd = c(); glm_diff_pca_sens = c(); glm_diff_pca_spec = c()
qda_diff_pca_accs = c(); qda_diff_pca_sd = c(); qda_diff_pca_sens = c(); qda_diff_pca_spec = c()
boo_diff_pca_accs = c(); boo_diff_pca_sd = c(); boo_diff_pca_sens = c(); boo_diff_pca_spec = c()

# Tests each amount of PCA dimensions.
for (diff in 2:8){
  
  lda_p4 = c(); lda_tru4 = c(); lda_acc4 = c()
  glm_p4 = c(); glm_tru4 = c(); glm_acc4 = c()
  qda_p4 = c(); qda_tru4 = c(); qda_acc4 = c()
  boo_p4 = c(); boo_tru4 = c(); boo_acc4 = c()
  
  # For 5-Fold Cross Validation.
  for (i in 1:k){
    samp = which(rand == i)
    
    if (COVARIATES == "nine"){
      train_mat = as.matrix(samp_data[-samp, 8:16])
      test_mat = as.matrix(samp_data[samp, 8:16])
    }else if (COVARIATES == "twelve"){
      if (VERSION == "MIT"){
        train_mat = as.matrix(samp_data[-samp, 8:19])
        test_mat = as.matrix(samp_data[samp, 8:19])
      }else{
        train_mat = as.matrix(samp_data[-samp, c(8:16, 18:20)])
        test_mat = as.matrix(samp_data[samp, c(8:16, 18:20)])
      }
    }
    
    cor_merged = prcomp(train_mat, retx = TRUE, center = TRUE, scale = TRUE)
    raw_merged = prcomp(train_mat, retx = TRUE, scale = FALSE) 
    
    train_pca = (train_mat) %*% raw_merged$rotation
    test_pca = test_mat %*% raw_merged$rotation
    
    cor_mat1 = cor(train_mat)
    lambdas = eigen(cor_mat1)$values
    var_exp = cumsum(lambdas) / sum(lambdas)
    
    PCA = diff
    
    train = as.data.frame(train_pca[, 1:PCA])
    
    test = as.data.frame(test_pca[, 1:PCA])
    
    train$State = samp_data[-samp, 5]
    test$State = samp_data[samp, 5]
    
    formula = State ~ .
    
    glm_res = glm_atr_funct(formula, train, test)
    lda_res = lda_atr_funct(formula, train, test)
    qda_res = qda_atr_funct(formula, train, test)
    boo_res = boo_atr_funct(formula, train, test)
    
    glm_acc4 = c(glm_acc4, glm_res[[1]]); glm_tru4 = c(glm_tru4, glm_res[[2]]); glm_p4 = c(glm_p4, glm_res[[3]])
    lda_acc4 = c(lda_acc4, lda_res[[1]]); lda_tru4 = c(lda_tru4, lda_res[[2]]); lda_p4 = c(lda_p4, lda_res[[3]])
    qda_acc4 = c(qda_acc4, qda_res[[1]]); qda_tru4 = c(qda_tru4, qda_res[[2]]); qda_p4 = c(qda_p4, qda_res[[3]])
    boo_acc4 = c(boo_acc4, boo_res[[1]]); boo_tru4 = c(boo_tru4, boo_res[[2]]); boo_p4 = c(boo_p4, boo_res[[3]])
  }
  
  glm_diff_pca_accs = c(glm_diff_pca_accs, mean(glm_acc4))
  glm_diff_pca_sd = c(glm_diff_pca_sd, sd(glm_acc4))
  glm_tab4 = table(glm_p4, glm_tru4); glm_tab4
  glm_diff_pca_sens = c(glm_diff_pca_sens, glm_tab4[1,1] / sum(glm_tab4[,1]))
  glm_diff_pca_spec = c(glm_diff_pca_spec, glm_tab4[2,2] / sum(glm_tab4[,2]))
  
  lda_diff_pca_accs = c(lda_diff_pca_accs, mean(lda_acc4))
  lda_diff_pca_sd = c(lda_diff_pca_sd, sd(glm_acc4))
  lda_tab4 = table(lda_p4, lda_tru4); lda_tab4
  lda_diff_pca_sens = c(lda_diff_pca_sens, lda_tab4[1,1] / sum(lda_tab4[,1]))
  lda_diff_pca_spec = c(lda_diff_pca_spec, lda_tab4[2,2] / sum(lda_tab4[,2]))
  
  qda_diff_pca_accs = c(qda_diff_pca_accs, mean(qda_acc4))
  qda_diff_pca_sd = c(qda_diff_pca_sd, sd(glm_acc4))
  qda_tab4 = table(qda_p4, qda_tru4); qda_tab4
  qda_diff_pca_sens = c(qda_diff_pca_sens, qda_tab4[1,1] / sum(qda_tab4[,1]))
  qda_diff_pca_spec = c(qda_diff_pca_spec, qda_tab4[2,2] / sum(qda_tab4[,2]))
  
  boo_diff_pca_accs = c(boo_diff_pca_accs, mean(boo_acc4))
  boo_diff_pca_sd = c(boo_diff_pca_sd, sd(glm_acc4))
  boo_tab4 = table(boo_p4, boo_tru4); boo_tab4
  boo_diff_pca_sens = c(boo_diff_pca_sens, boo_tab4[1,1] / sum(boo_tab4[,1]))
  boo_diff_pca_spec = c(boo_diff_pca_spec, boo_tab4[2,2] / sum(boo_tab4[,2]))
}
glm_diff_pca_accs; glm_diff_pca_sd; glm_diff_pca_sens; glm_diff_pca_spec
lda_diff_pca_accs; lda_diff_pca_sd; lda_diff_pca_sens; lda_diff_pca_spec
qda_diff_pca_accs; qda_diff_pca_sd; qda_diff_pca_sens; qda_diff_pca_spec
boo_diff_pca_accs; boo_diff_pca_sd; boo_diff_pca_sens; boo_diff_pca_spec

#### Results ####

## Model 1
nine_accs = c(); nine_sd = c(); nine_sens = c(); nine_spec = c()
for (i in 1:4){
  nine_accs = c(nine_accs, nine_cov_results[[i]][[1]])
  nine_sd = c(nine_sd, nine_cov_results[[i]][[2]])
  nine_sens = c(nine_sens, nine_cov_results[[i]][[4]])
  nine_spec = c(nine_spec, nine_cov_results[[i]][[5]])
}
paste("Nine-Transitions")
nine_accs; nine_sd; nine_sens; nine_spec

## Model 2
reg_long_accs = c(); reg_long_sd = c(); reg_long_sens = c(); reg_long_spec = c()
for (i in 1:4){
  reg_long_accs = c(reg_long_accs, reg_long_results[[i]][[1]])
  reg_long_sd = c(reg_long_sd, reg_long_results[[i]][[2]])
  reg_long_sens = c(reg_long_sens, reg_long_results[[i]][[4]])
  reg_long_spec = c(reg_long_spec, reg_long_results[[i]][[5]])
}
paste("Reg-to-Long")
reg_long_accs; reg_long_sd; reg_long_sens; reg_long_spec

## Model 3
reg_reg_accs = c(); reg_reg_sd = c(); reg_reg_sens = c(); reg_reg_spec = c()
for (i in 1:4){
  reg_reg_accs = c(reg_reg_accs, reg_reg_results[[i]][[1]])
  reg_reg_sd = c(reg_reg_sd, reg_reg_results[[i]][[2]])
  reg_reg_sens = c(reg_reg_sens, reg_reg_results[[i]][[4]])
  reg_reg_spec = c(reg_reg_spec, reg_reg_results[[i]][[5]])
}
paste("Reg-to-Reg")
reg_reg_accs; reg_reg_sd; reg_reg_sens; reg_reg_spec

## Model 4
pca_accs = c(max(glm_diff_pca_accs), max(lda_diff_pca_accs), max(qda_diff_pca_accs), max(boo_diff_pca_accs))
pca_indices = c(which(max(glm_diff_pca_accs) == glm_diff_pca_accs), which(max(lda_diff_pca_accs) == lda_diff_pca_accs), 
                which(max(qda_diff_pca_accs) == qda_diff_pca_accs), which(max(boo_diff_pca_accs) == boo_diff_pca_accs))
pca_sd = c(glm_diff_pca_sd[pca_indices[1]], lda_diff_pca_sd[pca_indices[2]], 
           qda_diff_pca_sd[pca_indices[3]], boo_diff_pca_sd[pca_indices[4]])
pca_sens = c(glm_diff_pca_sens[pca_indices[1]], lda_diff_pca_sens[pca_indices[2]], 
             qda_diff_pca_sens[pca_indices[3]], boo_diff_pca_sens[pca_indices[4]])
pca_spec = c(glm_diff_pca_spec[pca_indices[1]], lda_diff_pca_spec[pca_indices[2]], 
             qda_diff_pca_spec[pca_indices[3]], boo_diff_pca_spec[pca_indices[4]])
paste("PCA (Max)")
pca_accs; pca_sd; pca_sens; pca_spec; pca_indices

#### Graph Creation #### (Optional) ###############################################################

methods = c("Logistic", "LDA", "QDA", "Boosting")

## Model 1
nine_res = as.data.frame(matrix(c(methods, nine_accs, nine_sd, nine_sens, nine_spec), 4, 5, byrow = F))
colnames(nine_res) = c("Method", "Accuracy", "SE", "Sensitivity", "Specificity")
write.csv(nine_res, "C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/nine_results_1.csv")
nine_res = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/nine_results_1.csv", header = T)
nine_melt = melt(nine_res, id = "Method", measure = c("Accuracy", "SE", "Sensitivity", "Specificity"))

ggplot(nine_melt, aes(x = variable, y = value, fill = Method)) + geom_bar(stat = 'identity', position = 'dodge') + ylab("Value") + xlab("") +
  ggtitle("Response: AFib/Non-AFib | Covariates: Nine Transition States") + 
  geom_text(aes(label=round(value, 3)), position=position_dodge(0.9), size = 3.5, vjust = -0.5)

ggplot(nine_res, aes(x = Method, y = SE, fill = Method)) + geom_bar(stat = 'identity', position = 'dodge') + ylab("Standard Error") + xlab("") +
  ggtitle("Covariates: Nine Transition States") + 
  geom_text(aes(label=round(SE, 3)), position=position_dodge(0.9), size = 3.5, vjust = -0.5)

a = ggplot(nine_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Accuracy))
b = ggplot(nine_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = SE))
c = ggplot(nine_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Sens.))
d = ggplot(nine_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Spec.))
figure = ggarrange(a, b, c, d, labels = c(), ncol = 2, nrow = 2); figure

## Model 2
reg_long_res = as.data.frame(matrix(c(methods, reg_long_accs, reg_long_sd, reg_long_sens, reg_long_spec), 4, 5, byrow = F))
colnames(reg_long_res) = c("Method", "Accuracy", "SE", "Sensitivity", "Specificity")
write.csv(reg_long_res, "C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/reg_long_results_1.csv")
reg_long_res = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/reg_long_results_1.csv", header = T)
rl_melt = melt(reg_long_res, id = "Method", measure = c("Accuracy", "SE", "Sensitivity", "Specificity"))

ggplot(rl_melt, aes(x = variable, y = value, fill = Method)) + geom_bar(stat = 'identity', position = 'dodge') + ylab("Value") + xlab("") +
  ggtitle("Response: AFib/Non-AFib | Covariates: Regular-to-Long Transition State") + 
  geom_text(aes(label=round(value, 3)), position=position_dodge(0.9), size = 3.5, vjust = -0.5)

ggplot(reg_long_res, aes(x = Method, y = SE, fill = Method)) + geom_bar(stat = 'identity', position = 'dodge') + ylab("Standard Error") + xlab("") +
  ggtitle("Covariates: Regular-to-Long Transition State") + 
  geom_text(aes(label=round(SE, 3)), position=position_dodge(0.9), size = 3.5, vjust = -0.5)

a = ggplot(reg_long_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Accuracy))
b = ggplot(reg_long_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = SE))
c = ggplot(reg_long_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Sens.))
d = ggplot(reg_long_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Spec.))
figure = ggarrange(a, b, c, d, labels = c(), ncol = 2, nrow = 2); figure

## Model 3
reg_reg_res = as.data.frame(matrix(c(methods, reg_reg_accs, reg_reg_sd, reg_reg_sens, reg_reg_spec), 4, 5, byrow = F))
colnames(reg_reg_res) = c("Method", "Accuracy", "SE", "Sens.", "Spec.")
write.csv(reg_reg_res, "C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/reg_reg_results_1.csv")
reg_reg_res = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/reg_reg_results_1.csv", header = T)
a = ggplot(reg_reg_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Accuracy))
b = ggplot(reg_reg_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = SE))
c = ggplot(reg_reg_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Sens.))
d = ggplot(reg_reg_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Spec.))
figure = ggarrange(a, b, c, d, labels = c(), ncol = 2, nrow = 2); figure

## Model 4
pca_accs_res = as.data.frame(matrix(c(2:8, glm_diff_pca_accs, lda_diff_pca_accs, qda_diff_pca_accs, boo_diff_pca_accs), 7, 5, byrow = F))
colnames(pca_accs_res) = c("Dimensions", methods)
write.csv(pca_accs_res, "C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/pca_accs_results_1.csv")
pca_accs_res = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/pca_accs_results_1.csv", header = T)
a = ggplot(pca_accs_res, aes(x = Dimensions)) + geom_bar(stat = "identity", fill = hue_pal()(4)[3], aes(y = Logistic)) + guides(fill = F)
b = ggplot(pca_accs_res, aes(x = Dimensions)) + geom_bar(stat = "identity", fill = hue_pal()(4)[2], aes(y = LDA)) + guides(fill = F)
c = ggplot(pca_accs_res, aes(x = Dimensions)) + geom_bar(stat = "identity", fill = hue_pal()(4)[4], aes(y = QDA)) + guides(fill = F)
d = ggplot(pca_accs_res, aes(x = Dimensions)) + geom_bar(stat = "identity", fill = hue_pal()(4)[1], aes(y = Boosting)) + guides(fill = F)
figure = ggarrange(a, b, c, d, labels = c(), ncol = 2, nrow = 2); figure

pca_max_res = as.data.frame(matrix(c(methods, pca_accs, pca_indices + 1, pca_sd, pca_sens, pca_spec), 4, 6, byrow = F))
colnames(pca_max_res) = c("Method", "Accuracy", "Dim.", "SE", "Sens.", "Spec.")
write.csv(pca_max_res, "C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/pca_max_results_1.csv")

pca_max_res = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/pca_max_results_1.csv", header = T)
a = ggplot(pca_max_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Accuracy)) + 
  geom_text(aes(label=Dim., y = 0.05))
b = ggplot(pca_max_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = SE)) + 
  geom_text(aes(label=Dim., y = 0.0002))
c = ggplot(pca_max_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Sens.)) + 
  geom_text(aes(label=Dim., y = 0.05))
d = ggplot(pca_max_res, aes(x = Method)) + geom_bar(stat = "identity", aes(fill = Method, y = Spec.)) + 
  geom_text(aes(label=Dim., y = 0.05))
figure = ggarrange(a, b, c, d, labels = c(), ncol = 2, nrow = 2); figure

############################################################################################################################################