# Jericho Lawson
# Summer 2021
# Block Bootstrapping Testing

# libraries
library(boot); library(ggplot2); #library(plotly)

# MIT data
data = read.csv("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Data/all_trans_dataMIT.csv")
head(data)

# function for tsboot
test.fun = function(data) {
  return(cbind(data$Heartbeats, data$Reg.L))
}

# grapher
grapher = function(type, data, path, d_type){
  if (length(type) == 1){
    if (type == "N"){
      if (d_type == "pred"){
        ggplot(data, aes(x=count, y=Heartbeats, group = 1, color = pred)) + geom_line() +
          labs(color = "Type", x = "Segment", y = "Heartbeats") +
          scale_color_manual(labels = c("N"), values = c("dodgerblue1")) 
      } else {
        ggplot(data, aes(x=X, y=Heartbeats, group = 1, color = State)) + geom_line() +
          labs(color = "Type", x = "Segment", y = "Heartbeats") +
          scale_color_manual(labels = c("N"), values = c("dodgerblue1")) 
      }
    } else if (type == "AFIB"){
      if (d_type == "pred"){
        ggplot(data, aes(x=count, y=Heartbeats, group = 1, color = pred)) + geom_line() +
          labs(color = "Type", x = "Segment", y = "Heartbeats") +
          scale_color_manual(labels = c("AFIB"), values = c("tomato1")) 
      } else {
        ggplot(data, aes(x=X, y=Heartbeats, group = 1, color = State)) + geom_line() +
          labs(color = "Type", x = "Segment", y = "Heartbeats") +
          scale_color_manual(labels = c("AFIB"), values = c("tomato1")) 
      }
    }
  }
  else {
    if (d_type == "pred"){
      ggplot(data, aes(x=count, y=Heartbeats, group = 1, color = pred)) + geom_line() +
        labs(color = "Type", x = "Segment", y = "Heartbeats") +
        scale_color_manual(values = c("tomato1", "dodgerblue1")) 
    } else {
      ggplot(data, aes(x=X, y=Heartbeats, group = 1, color = State)) + geom_line() +
        labs(color = "Type", x = "Segment", y = "Heartbeats") +
        scale_color_manual(values = c("tomato1", "dodgerblue1")) 
    }
  }
  ggsave(path, width = 8, height = 4, units = "in")
}

subs = c("4015", "4043", "4048", "4126", 
         "4746", "4908", "4936", "5091", 
         "5121", "5261", "6426", "6453", 
         "6995", "7162", "7859", "7879", 
         "7910", "8215", "8219", "8378", 
         "8405", "8434", "8455")

bootstrap = function(sub, whole = F){
  if (whole == F){
    data_test = data[data$Subject == sub, ]
  } else{
    data_test = data
  }
  test = data_test[, c(6, 13)]
  data_test$X = data_test$X - data_test[1, ]$X + 1
  
  #### BLOCK BOOTSTRAPPING ####
  
  # create tallies for data
  i = 1; count = 0; tallies = data.frame()
  while(i <= nrow(data_test)){
    if (i == 1){
      curr = data_test$State[i]
      count = count + 1
    } else {
      if (data_test$State[i] != curr){
        row = c(curr, count)
        tallies = rbind(tallies, row)
        count = 1
        curr = data_test$State[i]
      } else {
        count = count + 1
      }
    }
    i = i + 1
  }
  row = c(curr, count)
  tallies = rbind(tallies, row)
  colnames(tallies) = c("State", "Tally")
  tallies$Tally = as.numeric(tallies$Tally)
  
  
  # the stationary bootstrap with mean block length 20
  set.seed(1)
  a = tsboot(test, test.fun, R = 100, 
             l = mean(tallies[tallies$Tally <= quantile(tallies$Tally, 0.75) + 1.5 * IQR(tallies$Tally), ]$Tally), 
             sim = "geom")
  
  test = data.frame("Heartbeats" = a$t[1, 1:dim(data_test)[1]])
  test$Reg.L = a$t[1, (dim(data_test)[1] + 1):(2 * dim(data_test)[1])]
  colnames(test) = c("Heartbeats", "Reg.L")
  
  # predictions for AFib/Non-AFib
  glm = glm(factor(State) ~ Reg.L, data = data_test, family = "binomial")
  glm_probs = predict(glm, test, type = "response")
  glm_pred = rep("AFIB", dim(test)[1]); glm_pred[glm_probs > .5] = "N"
  #list(mean(glm_pred == test$State), test$State, glm_pred)

  test$pred = glm_pred
  test$count = 1:length(test$pred)

  sink(file = paste("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Graphs/Bootstrapping/sub_", sub, "comp.txt", sep =""))
  print("Tallies for Pred.")
  print(table(test$pred))
  print("Tallies for Data")
  print(table(data_test$State))
  sink(file = NULL)

  path = paste("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Graphs/Bootstrapping/sub_", sub, "block.jpg", sep ="")
  grapher(levels(factor(glm_pred)), test, path, "pred")

  path = paste("C:/Users/Mario/Documents/Internship and REU Stuff/UNCW 2019/Atrial_Fibrillation/Summer 2021 Testing/Graphs/Bootstrapping/sub_", sub, "real.jpg", sep ="")
  grapher(levels(factor(data_test$State)), data_test, path, "test")
}

sapply(subs, bootstrap)

sapply(c("whole"), bootstrap, whole = T)


##############

# plots to compare
par(mfrow = c(2, 1))
plot(a$t[1, ], type = "l", main = 'bootstrap1')
plot(data$Reg.L, type = "l", main = 'standard')

ggplot(data_test, aes(x=X, y=Reg.L, group = 1, color = State)) + geom_line()

# characteristics of tallies
hist(tallies$Tally, breaks = 50)
plot(density(tallies[tallies$Tally <= 59.5, ]$Tally))
lines(dgeom(0:60, prob = 1 / mean(tallies[tallies$Tally <= 59.5, ]$Tally)), col = "red")

mean(tallies$Tally)

quantile(tallies$Tally, 0.75) + 1.5 * IQR(tallies$Tally)

boxplot(tallies$Tally)

mean(tallies[tallies$Tally <= 59.5, ]$Tally)
