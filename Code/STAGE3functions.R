## Jericho Lawson
## Summer 2019, 2021
## AFib Functions for Training Program

# Functions that aid the STAGE3training.R code to process the data, train
# various models, and output success metrics.

#### Utility Functions #########################################################

# Function to open a specific library depending on which algorithms are used.
# I: pick (3-letter code indicating the algorithm in use)
# O: none
library_selector = function(pick){
  if (pick == "LDA" | pick == "QDA"){ # LDA/QDA
    library(MASS)
  } else if (pick == "GBM"){ # gradient boosting
    library(gbm)
  } else if (pick == "XGB"){ # XGBoost
    library(xgboost)
  } else if (pick == "LGB"){ # light graident boosting (LightGBM)
    library(lightgbm)
  } else if (pick == "SVM"){ # support-vector machine
    library(e1071)
  } else if (pick == "RFO"){ # random forest
    library(randomForest)
  }
}

# Creates formula based on user-specified model.
# I: indices (indices of covariates to be used from data frame)
#    cols (names of columns from data frame)
# O: formula type that has the model that will be used
formula_maker = function(indices, cols){
  names = cols[indices]
  return(formula(paste("State ~", paste(names, collapse = "+"))))
}

# Labels subjects (randomly for 2017/CinC dataset, by subject for MIT-BIH)
# for k-fold cross validation.
# I: data (data frame in use), vers ("M" or "C"),
#    seed (user-generated seed for randomization), 
#    k (number of folds, 2017 data only)
#    cols (names of columns from data frame)
# O: Numeric vector that labels which observation is in which fold for
#    cross-validation purposes
grouper = function(data, vers, seed, k){
  if (vers == "M"){ # MIT-BIH data
    return(as.numeric(factor(data$Subject)))
  } else { # 2017/CinC data
    set.seed(seed)
    return(sample(rep(1:k, length.out = dim(data)[1])))
  }
}

# Chooses which function to run based on the user-specified algorithm.
# I: pick (3-letter code indicating the algorithm in use),
#    formula (model used for algorithm),
#    index_c (indices in data frame that correspond to covariate data),
#    index_r (index in data frame that corresponds to response data),
#    train (training data),
#    test (testing data),
#    thresh (probabilistic cutoff point that classifies segment as AFib or N),
#    seed (seed for random generation)
# O: Results from algorithm (list containing true/predicted states and time)
algorithm_picker = function(pick, formula, index_c, index_r, train, 
                            test, thresh, seed){
  if (pick == "LOG"){ # logistic regression
    logistic(formula, train, test, thresh, seed)
  } else if (pick == "LDA"){ # LDA
    lda_f(formula, train, test, seed)
  } else if (pick == "QDA"){ # QDA
    qda_f(formula, train, test, seed)
  } else if (pick == "GBM"){ # gradient boosting
    boost_gbm(formula, train, test, thresh, seed)
  } else if (pick == "XGB"){ # XGBoost
    boost_xgb(index_c, index_r, train, test, 
              nrounds = 5000, eta = 0.01, gamma = 0, 
              max_depth = 4, min_child_weight = 1, subsample = 1, 
              colsample_bytree = 1, lambda = 0, alpha = 1, 
              thresh = 0.5, seed = seed)
  } else if (pick == "LGB"){ # light gradient boosting (LightGBM)
    boost_lgb(index_c, index_r, train, test, nrounds = 100, type = "gbdt", 
              num_leaves = 31, max_depth = -1, 
              num_threads = 0, thresh, seed)
  } else if (pick == "SVM"){ # support-vector machine
    svm_f(formula, index_c, train, test, kernel = "linear", cost = 1, seed)
  } else if (pick == "RFO"){ # random forest
    rf_f(formula, train, test, ntrees = 500, seed)
  }
}

# Runs all testing for k-fold cross validation and times the process. 
# Returns true/predicted states of each fold, as well as how long it took 
# to run all folds combined (in seconds).
# I: pick (3-letter code indicating the algorithm in use),
#    data (data set used)
#    rand (vector of labels specifying fold for cross-validation purposes)
#    formula (model used for algorithm),
#    index_c (indices in data frame that correspond to covariate data),
#    index_r (index in data frame that corresponds to response data),
#    thresh (probabilistic cutoff point that classifies segment as AFib or N),
#    seed (seed for random generation),
#    folds (number of folds needed for k-fold cross validation)
# O: List containing true/predicted states of observation in each fold in
#    the first element and the time it took to run k-fold cross validation
#    in the second element.
testing = function(pick, data, rand, formula, index_c, index_r, thresh, 
                   seed, folds = max(rand)){
  # create empty list for results and records beginning time
  results = list()
  time = Sys.time()
  
  for (i in 1:folds){ # for each fold
    # creates training/testing sets
    samp = which(rand == i); train = data[-samp, ]; test = data[samp, ] 
    
    # runs algorithm, assigns results to ith element of list
    results[[i]] = algorithm_picker(pick, formula = formula, index_c = index_c, 
                                    index_r = index_r, train = train, 
                                    test = test, thresh = 0.5, seed = seed)
  }
  
  # returns list of results and the length of time involved in running this 
  # process
  return(list("results" = results, 
              "time" = as.numeric(difftime(Sys.time(), time, units = "secs"))))
}

# Using the list of results and times, the function outputs various metrics
# for verifying the success of a model, such as sensitivity, prediction
# accuracy, and F-1 score.
# I: list (contains true/predicted states of observation in each fold in
#          the first element and the time it took to run k-fold cross 
#          validation in the second element)
# O: List containing true/predicted states of observation in each fold in
#    the first element and the time it took to run k-fold cross validation
#    in the second element.
metrics_f = function(list){
  # create confusion matrix
  for (fold in 1:length(list$results)){
    if (fold == 1){
      confusion = table(list$results[[fold]]$Truth, list$results[[fold]]$Pred)
    } else {
      confusion = confusion + table(list$results[[fold]]$Truth, list$results[[fold]]$Pred)
    }
  }
  
  # Returns list of metrics.
  # Note: works when there are tallies of zero!
  return(list("confusion" = confusion, 
              "sensitivity" = confusion[2, 2] / sum(confusion[2, ]),
              "specificity" = confusion[1, 1] / sum(confusion[1, ]),
              "f1" = confusion[2, 2] / (confusion[2, 2] + 0.5 * (confusion[1, 2] + confusion[2, 1])),
              "accuracy" = (confusion[1, 1] + confusion[2, 2]) / sum(confusion),
              "time" = list$time))
}

#### Algorithms ################################################################

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
  return(list("Truth" = factor(test$State, levels = c("N", "AFIB")), "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# LightGBM testing.
boost_lgb = function(index_c, index_r, train, test, nrounds = 100, type = "gbdt", 
                     num_leaves = 31, max_depth = -1, num_threads = 0, thresh, seed){
  # Reclassifies states such that they are binary.
  train$State = ifelse(train$State == "AFIB", 1, 0)
  test$State = ifelse(test$State == "AFIB", 1, 0)
  
  # Converts training set to a lgb.Dataset object for training.
  dtrain = lgb.Dataset(data = as.matrix(train[, index_c]), label = train[, index_r])
  
  # Trains model with training data.
  set.seed(seed)
  mod = lgb.train(data = dtrain, boosting = type, num_leaves = num_leaves, nrounds = nrounds, 
                  objective = "binary", max_depth = max_depth, verbose = 0, num_threads = num_threads,
                  force_row_wise = T)
  pred = predict(mod, as.matrix(test[, index_c]))
  pred = ifelse(pred >= thresh, "AFIB", "N")
  
  test$State = ifelse(test[, index_r] == 1, "AFIB", "N")
  
  # Returns actual and predicted states.
  return(list("Truth" = factor(test$State, levels = c("N", "AFIB")), "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# XGBoost testing.
boost_xgb = function(index_c, index_r, train, test, nrounds, eta, gamma, max_depth, min_child_weight, subsample, colsample_bytree, lambda, alpha, thresh, seed){
  # Reclassifies states such that they are binary. Converts data frames to matrices.
  train$State = ifelse(train$State == "AFIB", 1, 0); train = data.matrix(train)
  test$State = ifelse(test$State == "AFIB", 1, 0); test = data.matrix(test)
  
  # Trains model with training data.
  set.seed(seed)
  mod = xgboost(data = as.matrix(train[, index_c]), label = train[, index_r], 
                nrounds = nrounds, eta = eta, gamma = gamma,
                max_depth = max_depth, min_child_weight = min_child_weight, 
                subsample = subsample, colsample_bytree = colsample_bytree, 
                lambda = lambda, alpha = alpha, 
                scale_pos_weight = sum(train[, index_r] == 0) / sum(train[, index_r] == 1),
                objective = "reg:squarederror", verbose = 0)
  pred = predict(mod, as.matrix(test[, index_c]))
  pred = ifelse(pred >= thresh, "AFIB", "N")
  
  test[, index_r] = ifelse(test[, index_r] == 1, "AFIB", "N")
  
  # Returns actual and predicted states.
  return(list("Truth" = factor(test[, index_r], levels = c("N", "AFIB")), "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Support-Vector Machine testing.
svm_f = function(formula, index_c, train, test, kernel = "radial", cost = 1, seed){
  # Trains model with training data.
  set.seed(seed)
  mod = svm(form, data = train, kernel = kernel, cost = cost, scale = FALSE)
  
  # Classifies test data based on training model.
  pred = predict(mod, test[, index_c])

  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = factor(pred, levels = c("N", "AFIB"))))
}

# Random forest testing.
# Linear discriminant analysis testing.
rf_f = function(formula, train, test, ntrees, seed){
  # Trains model with training data.
  set.seed(seed)
  mod = randomForest(form, data = train, ntree = ntrees)
  
  # Classifies test data based on training model.
  pred = predict(mod, test)
  
  # Returns actual and predicted states.
  return(list("Truth" = test$State, "Pred" = pred))
}

################################################################################