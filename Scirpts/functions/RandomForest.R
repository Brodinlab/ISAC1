# Functions to perform random forest analysis

# data: scaled data matrix
# group: response vector corresponding to the data (reference group as 0, case group as 1)

# Random forest step 1: model fit and variable selection -----

rf_fit <- function(data, group){
  results <- list()
  
  # Model fit and prediction
  print("Fitting model...")
  
  # Grid of parameter combinations
  tgrid <- expand.grid(
    .mtry = seq(from = 5, to = 200, by = 40), # number of variables randomly selected by each tree
    .splitrule = "gini", # split rule
    .min.node.size = seq(from = 2, to = 32, by = 5) # min node size
  )
  
  # Convert group as factor
  groupfact  <- as.factor(ifelse(group == 1,"B","A"))
  my.data <- data.frame(groupfact, data) 
  
  # Train the model
  set.seed(1109)
  RF.fit <- train(groupfact  ~ ., data = my.data,
                  method = "ranger", # random forest
                  trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE, classProbs = TRUE, savePredictions = TRUE), # 5-fold cross validation
                  tuneGrid = tgrid, # parameter grid
                  num.trees = 100, # number of trees in the random forest
                  importance = "permutation") # use permutation to calculate variable importance
  
  # Prediction using trained model
  RF.pred <- predict(RF.fit, data)
  
  # Confusion matrix
  print("Prediction results of fitted model:")
  print(table(groupfact, RF.pred))
  
  # Calculate and plot ROC curve
  selected_indices <- RF.fit$pred$mtry == RF.fit$bestTune$mtry
  print("ROC curve:")
  roc_curve <- plot.roc(RF.fit$pred$obs[selected_indices], RF.fit$pred$B[selected_indices], print.auc = TRUE)
  
  print(paste("Repeated fitting and variable selection..."))
  # Create a matrix to store results of variable importance
  RF.imp.mat <- matrix(0, ncol=100, nrow=ncol(data))
  # Number of repetition
  n.repeated <- 10
  # Repeated fitting
  set.seed(1109)
  for(i in 1:n.repeated){
    RF.fit <- train(groupfact  ~ ., data = my.data,
                    method = "ranger",
                    trControl = trainControl(method="cv", number = 5, verboseIter = FALSE, classProbs = FALSE),
                    tuneGrid = tgrid,
                    num.trees = 100,
                    importance = "permutation")
    RF.imp.mat[,i] <- (varImp(RF.fit))$importance$Overall # save variable importance of each fitting
  }
  
  # Calculate mean importance of each variable in all repetitions
  RF.imp <- apply(RF.imp.mat, 1, function(x) {mean(x)})
  RF.selected.proteins <- which(RF.imp > sort(RF.imp,decreasing = TRUE)[31])  # top 30 proteins
  print("Selected features:")
  print(colnames(data)[RF.selected.proteins])
  
  # Save results
  results$index.selected.proteins <- RF.selected.proteins
  results$names.selected.proteins <- colnames(data)[RF.selected.proteins]
  results$importances <- RF.imp
  results$tgrid <- tgrid
  results$roc_curve <- roc_curve
  
  return(results)
}

# Random forest step 2: evaluating performance (nested CV) -----

rf_evaluate <- function(data, group, res_rf){ # results from rf_fit
  results <- res_rf
  print("Evaluating performance...")
  
  # Performance: CV
  set.seed(1109)
  ncv <- 100
  cv.sens <- cv.spec <- cv.ppv <- cv.npv <- cv.auc <- c()
  TP <- TN <- FP <- FN <- c()
  
  for(cv in 1:ncv){ # outer CV loop
    if(cv %% 10 == 0) print(cv)
    set.seed(1109)
    cv.index <- createDataPartition(y = as.factor(group), p = 0.1, list = FALSE) # 100 outer CVs with 10% test sets, stratified! as using y as factor
    x.train <- data[-cv.index,]
    y.train <- group[-cv.index]
    x.test  <- data[cv.index,]
    y.test  <- ifelse(group[cv.index] == 1, "B", "A")
    
    groupfact <- as.factor(ifelse(y.train == 1, "B", "A"))
    my.data <- data.frame(groupfact, x.train) 
    set.seed(1109)
    RF.fit <- train(groupfact  ~ ., data = my.data,
                    method = "ranger",
                    trControl = trainControl(method="cv", number = 5, verboseIter = FALSE, classProbs = TRUE),
                    tuneGrid = results$tgrid,
                    num.trees = 100,
                    importance = "permutation")
    
    pred.prob <- (predict(RF.fit, x.test, type="prob"))$B
    pred <- predict(RF.fit, x.test)
    
    test.res <- cbind(c(sum(y.test == "B" & pred == "B"),sum(y.test == "A" & pred == "B")),
                      c(sum(y.test == "B" & pred == "A"),sum(y.test == "A" & pred == "A")))
    
    TP[cv] = test.res[1,1]; TN[cv] = test.res[2,2]; FP[cv] = test.res[2,1]; FN[cv] = test.res[1,2]
    
    cv.sens[cv] <- TP[cv]/(TP[cv]+FN[cv])
    cv.spec[cv] <- TN[cv]/(TN[cv]+FP[cv])
    cv.ppv[cv]  <- TP[cv]/(TP[cv]+FP[cv])
    cv.npv[cv]  <- TN[cv]/(TN[cv]+FN[cv])
    cv.auc[cv]  <-  as.numeric(pROC::auc(response = y.test, predictor = pred.prob, quiet = TRUE))
  }
  
  # Figures for performance results:
  print("Performance results of fitted model: plot")
  par(mfrow = c(2,2))
  try({
    hist(cv.sens,main=paste("sens",paste(round(quantile(cv.sens,probs = c(0.025,0.5,0.975), na.rm = TRUE),digits = 4),collapse=" ")))
    hist(cv.spec,main=paste("spec",paste(round(quantile(cv.spec,probs = c(0.025,0.5,0.975), na.rm = TRUE),digits = 4),collapse=" ")))
    hist(cv.ppv, main=paste("ppv", paste(round(quantile(cv.ppv,probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 4),collapse=" ")))
    hist(cv.npv, main=paste("npv", paste(round(quantile(cv.npv,probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 4),collapse=" ")))
  }, silent = TRUE)
  par(mfrow = c(1,1))
  
  # Mean confusion matrix:
  tab <- cbind(c(mean(TP), mean(FP)), c(mean(FN), mean(TN)))
  rownames(tab) <- c("1","0")
  colnames(tab) <- c("pred 1","pred 0")
  print("Mean confusion matrix:")
  print(tab)
  
  results$TP <- TP
  results$TN <- TN
  results$FP <- FP
  results$FN <- FN
  results$cv.sens <- cv.sens
  results$cv.spec <- cv.spec
  results$cv.ppv <- cv.ppv
  results$cv.npv <- cv.npv
  results$cv.auc <- cv.auc
  
  return(results)
}

# Random forest step 3: establishing significance of original performance (nested CV) -----

rf_significance <- function(data, group, res_rf){ # results from rf_evaluate
  results <- res_rf
  # Permuted labels performance (CV)
  print(paste("Permutation tests..."))
  
  set.seed(1109)
  ncv <- 100
  cv.sens.perm <- cv.spec.perm <- cv.ppv.perm <- cv.npv.perm <- cv.auc.perm <- c()
  TP.perm <- TN.perm <- FP.perm <- FN.perm <- c()
  
  for(cv in 1:ncv){ # outer CV loop 
    group.perm <- sample(group) 
    if(cv %% 10 == 0) print(cv)
    set.seed(1109)
    cv.index <- createDataPartition(y = as.factor(group.perm), p = 0.1, list=FALSE) # stratified partitions
    x.train <- data[-cv.index,]
    y.train <- group.perm[-cv.index]
    x.test  <- data[cv.index,]
    y.test  <- ifelse(group.perm[cv.index] == 1, "B", "A")
    
    groupfact  <- as.factor(ifelse(y.train == 1, "B", "A"))
    my.data <- data.frame(groupfact, x.train) 
    set.seed(1109)
    RF.fit <- train(groupfact  ~ ., data = my.data,
                    method = "ranger",
                    trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE, classProbs = TRUE),
                    tuneGrid = results$tgrid,
                    num.trees = 100,
                    importance = "permutation")
    
    pred.prob <- (predict(RF.fit, x.test, type="prob"))$B
    pred <- predict(RF.fit, x.test)
    
    test.res <- cbind(c(sum(y.test == "B" & pred == "B"),sum(y.test == "A" & pred == "B")),
                      c(sum(y.test == "B" & pred == "A"),sum(y.test == "A" & pred == "A")))
    
    TP.perm[cv] = test.res[1,1]; TN.perm[cv] = test.res[2,2]; FP.perm[cv] = test.res[2,1]; FN.perm[cv] = test.res[1,2]
    
    cv.sens.perm[cv] <- TP.perm[cv]/(TP.perm[cv]+FN.perm[cv])
    cv.spec.perm[cv] <- TN.perm[cv]/(TN.perm[cv]+FP.perm[cv])
    cv.ppv.perm[cv]  <- TP.perm[cv]/(TP.perm[cv]+FP.perm[cv])
    cv.npv.perm[cv]  <- TN.perm[cv]/(TN.perm[cv]+FN.perm[cv])
    cv.auc.perm[cv]  <- as.numeric(pROC::auc(response=y.test, predictor = pred.prob, quiet=TRUE))
  }
  
  # Figures for performance results:
  print("Performance results of permuation test: plot")
  par(mfrow = c(2,2))
  try({
    hist(cv.sens.perm,main = paste("perm sens",paste(round(quantile(cv.sens.perm, probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 2), collapse = " ")))
    hist(cv.spec.perm,main = paste("perm spec",paste(round(quantile(cv.spec.perm, probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 2), collapse = " ")))
    hist(cv.ppv.perm,main = paste("perm ppv",paste(round(quantile(cv.ppv.perm, probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 2), collapse = " ")))
    hist(cv.npv.perm,main = paste("perm npv",paste(round(quantile(cv.npv.perm, probs = c(0.025,0.5,0.975), na.rm = TRUE), digits = 2), collapse = " ")))
  }, silent = TRUE)
  par(mfrow = c(1,1))
  
  # Mean confusion matrix (permutation test):
  tab <- cbind(c(mean(TP.perm), mean(FP.perm)), c(mean(FN.perm), mean(TN.perm)))
  rownames(tab) <- c("0", "1")
  colnames(tab) <- c("pred 1", "pred 0")
  print("Mean confusion matrix of permutation test:")
  print(tab)
  
  results$TP.perm <- TP.perm
  results$TN.perm <- TN.perm
  results$FP.perm <- FP.perm
  results$FN.perm <- FN.perm
  results$cv.sens.perm <- cv.sens.perm
  results$cv.spec.perm <- cv.spec.perm
  results$cv.ppv.perm <- cv.ppv.perm
  results$cv.npv.perm <- cv.npv.perm
  results$cv.auc.perm <- cv.auc.perm
  
  results$p.val.perm.sens <- sum(cv.sens.perm > results$cv.sens,na.rm=TRUE) / length(results$cv.sens)
  results$p.val.perm.spec <- sum(cv.spec.perm > results$cv.sens,na.rm=TRUE) / length(results$cv.spec)
  results$p.val.perm.ppv  <- sum(cv.ppv.perm  > results$cv.ppv,na.rm=TRUE) / length(results$cv.ppv)
  results$p.val.perm.npv  <- sum(cv.npv.perm  > results$cv.npv,na.rm=TRUE) / length(results$cv.npv)
  results$p.val.perm.auc  <- sum(cv.auc.perm  > results$cv.auc,na.rm=TRUE) / length(results$cv.auc)
  
  print(paste("p-value of sensitivity:", results$p.val.perm.sens))
  print(paste("p-value of specificity:", results$p.val.perm.spec))
  print(paste("p-value of positive prediction value:", results$p.val.perm.ppv))
  print(paste("p-value of negative prediction value:", results$p.val.perm.npv))
  print(paste("p-value of AUC:", results$p.val.perm.auc))
  
  return(results)
}

# Plot random forest -----

plotRF <- function(res.rf, res.test.all){
  imp.df <- data.frame(feature = res.rf$names.selected.proteins, imp = res.rf$importances[res.rf$index.selected.proteins]) %>% 
    left_join(res.test.all %>% select(feature, log2fc, type)) %>%
    mutate(feature = fct_reorder(feature, imp)) %>% 
    mutate(direction = ifelse(log2fc > 0, "up", "down")) %>%
    arrange(desc(feature))
  
  p <- ggplot(imp.df, aes(x = feature, y = imp)) +
    geom_segment(aes(x = feature, xend = feature, y = 0, yend = imp, color = direction), alpha = 0.5, show.legend = FALSE) +
    geom_point(aes(shape = type, color = direction, fill = direction), size=4, alpha=0.7) +
    scale_color_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
    scale_fill_manual(values = c("down" = "grey30", "up" = "#BC80BDFF")) +
    scale_shape_manual(values = c("cell" = 21, "protein" = 23, "age" = 22)) +
    theme_light() +
    coord_flip() +
    labs(x = "Feature", y = "Permutation importance", color = "Direction", shape = "Type", title = title) +
    guides(fill = "none") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_line(linetype = "dashed"), 
      panel.grid.major = element_line(linetype = "dashed"), 
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = rev(c("down" = "black", "up" = "#BC80BDFF")[imp.df$direction])),
      plot.title = element_text(size = 10)
    )
  
  return(p)
}