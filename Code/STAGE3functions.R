## Jericho Lawson
## Summer 2019, 2021
## AFib Functions for Training Program



#### LIBRARIES #################################################################

packages = c("MASS", "gbm", "scales", "ggplot2", "reshape2", "ggpubr", "tidyr")
lapply(packages, library, character.only = TRUE)

formula_maker = function(indices, cols){
  names = cols[indices]
  return(formula(paste("State ~", paste(names, collapse = "+"))))
}

grouper = function(data, vers, seed, k){
  if (vers == "M"){
    return(as.numeric(factor(data$Subject)))
  } else {
    set.seed(seed)
    return(sample(rep(1:k, length.out = dim(data)[1])))
  }
}

#"LOG", "LDA", "QDA", "GBM"
algorithm_picker = function(pick, formula, train, test, thresh, seed){
  if (pick == "LOG"){
    logistic(formula, train, test, thresh, seed)
  } else if (pick == "LDA"){
    lda_f(formula, train, test, seed)
  } else if (pick == "QDA"){
    qda_f(formula, train, test, seed)
  } else if (pick == "GBM"){
    boost_gbm(formula, train, test, thresh, seed)
  } 
}

# Function that runs all testing.
testing = function(pick, data, rand, formula, thresh, seed, folds = max(rand)){
  results = list()
  
  for (i in 1:folds){
    samp = which(rand == i); train = data[-samp, ]; test = data[samp, ]
    results[[i]] = algorithm_picker(pick, formula = formula, train = train, test = test, thresh = 0.5, seed = seed)
  }
  
  return(results)
}

metrics = function(list){
  for (fold in 1:length(list)){
    if (fold == 1){
      confusion = table(list[[fold]]$Truth, list[[fold]]$Pred)
    } else {
      confusion = confusion + table(list[[fold]]$Truth, list[[fold]]$Pred)
    }
  }
  
  # works when there are tallies of zero!
  return(list("confusion" = confusion, 
              "sensitivity" = confusion[2, 2] / sum(confusion[2, ]),
              "specificity" = confusion[1, 1] / sum(confusion[1, ]),
              "f1" = confusion[2, 2] / (confusion[2, 2] + 0.5 * (confusion[1, 2] + confusion[2, 1])),
              "accuracy" = (confusion[1, 1] + confusion[2, 2]) / sum(confusion)))
}

#### ALGORITHMS ################################################################

# Logistic regression testing.
logistic = function(formula, train, test, thresh, seed){
  # Trains model with training data.
  set.seed(seed)
  mod = glm(formula, data = train, family = "binomial")
  
  # Classifies test data based on training model.
  probs = predict(mod, test, type = "response")
  pred = rep("N", dim(test)[1])
  pred[probs > thresh] = "AFIB"
  
  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Linear discriminant analysis testing.
lda_f = function(formula, train, test, seed){
  # Trains model with training data.
  set.seed(seed)
  mod = lda(formula, data = train)
  
  # Classifies test data based on training model.
  pred = predict(mod, test, type = "response")$class
  
  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Quadratic discriminant analysis testing.
qda_f = function(formula, train, test, seed){
  # Trains model with training data.
  set.seed(seed)
  mod = qda(formula, data = train)
  
  # Classifies test data based on training model.
  pred = predict(mod, test, type = "response")$class
  
  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Boosting testing.
boost_gbm = function(formula, train, test, thresh, seed, ntrees = 5000, shrink = 0.01, in_depth = 4){
  # Reclassifies states such that they are binary.
  train$State = ifelse(train$State == "AFIB", 1, 0)
  test$State = ifelse(test$State == "AFIB", 1, 0)
  
  # Trains model with training data.
  set.seed(seed)
  mod = gbm(formula, data = train, distribution = "bernoulli", n.trees = ntrees, shrinkage = shrink, interaction.depth = in_depth)
  
  # Classifies test data based on training model.
  pred = predict(mod, test, n.trees = 5000, type = "response")
  pred = ifelse(pred >= thresh, "AFIB", "N")
  
  # Reclassifies states back to AFib/N.
  train$State = ifelse(train$State == 1, "AFIB", "N")
  test$State = ifelse(test$State == 1, "AFIB", "N")
  
  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Currently working on items below ###########################

# XGBoost testing.
boost_xgb = function(index_r, index_c, train, test, nrounds, eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, lambda, alpha){
  # Reclassifies states such that they are binary. Converts data frames to matrices.
  train$State = ifelse(train$State == "AFIB", 1, 0); train = data.matrix(train)
  test$State = ifelse(test$State == "AFIB", 1, 0); test = data.matrix(test)
  
  # Trains model with training data.
  set.seed(seed)
  xgb = xgboost(data = as.matrix(train[, index_c]), label = train[, 5], 
                nrounds = nrounds, eta = eta, gamma = gamma,
                max_depth = max_depth, min_child_weight = min_child_weight, 
                subsample = subsample, colsample_bytree = colsample_bytree, 
                lambda = lambda, alpha = alpha, 
                scale_pos_weight = sum(train[, 5] == 0) / sum(train[, 5] == 1),
                objective = "binary:logistic", verbose = 0)
  xgb_pred = predict(xgb, as.matrix(test[, index]))
  xgb_pred = ifelse(xgb_pred >= 0.5, "AFIB", "N")
  
  test[,5] = ifelse(test[,5] == 1, "AFIB", "N")
  
  return(list(mean(xgb_pred == test[,5]), test[,5], xgb_pred, Sys.time() - start))
}

