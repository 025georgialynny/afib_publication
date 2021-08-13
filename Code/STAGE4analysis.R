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