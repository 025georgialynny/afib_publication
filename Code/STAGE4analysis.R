## Jericho Lawson
## Summer 2019, 2021
## Analysis for Results from Stage 3 Code

library(ggplot2)
library(gridExtra)

DIR = "C:/Users/mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/" # directory
setwd(DIR)

FILE = "Results/results_alltimes_8to16_all8_MIT.txt"

text = read.delim(FILE, header = F)
codes = c("LOG", "LDA", "QDA", "GBM", "XGB", "LGB", "SVM", "RFO")

results = list()
line = 1
ind = 0
extract = F
dataframe = F

info1 = c()
info2 = c()
info3 = c()
info4 = c()

while (line <= dim(text)[1]){

  if (substr(text[line, 1], 1, 18) == "[1] Segment Length"){
    
    if (length(results) == 0){
      extract = T
    }
    
    ind = ind + 1
    results[[ind]] = list()
    
    info1 = c(info1, rep(substr(text[line, 1], 21, nchar(text[line, 1])), 8))
    info2 = c(info2, rep(substr(text[line + 1, 1], 17, nchar(text[line + 1, 1])), 8))
    info3 = c(info3, rep(substr(text[line + 2, 1], 14, nchar(text[line + 2, 1])), 8))
    
    line = line + 3
  } else if (extract == T & nchar(text[line, 1]) == 4 & substr(text[line, 1], 2, 4) %in% codes){
    
    info4 = c(info4, substr(text[line, 1], 2, 4))
    
    metrics = c()
    
    # will add mechanism to extract confusion matrix at later time
    
    line = line + 6
    
    for (add in seq(0, 8, 2)){
      metrics = c(metrics, as.numeric(substr(text[line + add + 1, 1], 5, nchar(text[line + add + 1, 1]))))
    }
    
    line = line + 10
    
    if (dataframe == F){
      data = as.data.frame(matrix(metrics, nrow = 1, ncol = length(metrics)))
      dataframe = T
    } else {
      data = rbind(data, metrics)
    }

  } else {
    line = line + 1
  }
}

data = cbind(info4, info1, info2, info3, data)
colnames(data) = c("algorithm", "length", "covariates", "datatype", "sensitivity", "specificity", "f1", "accuracy", "comptime")


# comptime comparison
p1 = ggplot(data, aes(x = length, y = comptime, group = algorithm)) + geom_line(aes(color = algorithm)) + geom_point(aes(color = algorithm))

# sensitivity comparison
p2 = ggplot(data, aes(x = length, y = sensitivity, group = algorithm)) + geom_line(aes(color = algorithm)) + geom_point(aes(color = algorithm))

# specificity comparison
p3 = ggplot(data, aes(x = length, y = specificity, group = algorithm)) + geom_line(aes(color = algorithm)) + geom_point(aes(color = algorithm))

# f1 score comparison
p4 = ggplot(data, aes(x = length, y = f1, group = algorithm)) + geom_line(aes(color = algorithm)) + geom_point(aes(color = algorithm))

# accuracy comparison
p5 = ggplot(data, aes(x = length, y = accuracy, group = algorithm)) + geom_line(aes(color = algorithm)) + geom_point(aes(color = algorithm))

grid.arrange(p1, p2, p3, p4, p5, nrow = 3)




p1
p2
p3
p4
p5


############################ DO NOT RUN CODE BELOW ##################################################################################################


# originally taken from STAGE3training.R

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