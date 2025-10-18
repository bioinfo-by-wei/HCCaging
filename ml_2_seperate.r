rm(list = ls())
library(qs)
library("e1071")
library(pROC)
library(survcomp)
library(survival)
library(dplyr)
library(tidymodels)
library(xgboost)
library(randomForest)
library(gbm)
library(ggplot2)

plot_survival_and_timeROC <- function(train_data, pred_prob, store_name) {
  # Grouping: Divide risk values into high and low groups based on median
  pred_prob <- as.numeric(pred_prob)
  risk_score_plot <- ifelse(pred_prob > median(pred_prob), "High", "Low")
  # Construct data frame
  risk_data_plot <- data.frame(
    time = train_data$time,
    event = train_data$event,
    risk_score_plot = factor(
      risk_score_plot,
      levels = c("High", "Low"),
      labels = c("High", "Low")
    )
  )

  # Survival analysis model
  fit <- survfit(Surv(time, as.numeric(event)) ~ risk_score_plot, data = risk_data_plot)

  # Plot KM curve
  lasso_km <- ggsurvplot(
    fit,
    data = risk_data_plot,
    legend.labs = c("High", "Low"),
    pval = TRUE,
    pval.size = 8,
    risk.table = TRUE,
    size = 0.6,  # Line thickness
    risk.table.height = 0.25, # Ensure plotting on same page
    surv.median.line = "hv",
    legend.title = "RiskScore",
    title = "Overall survival",
    ylab = "Survival probability",
    xlab = "Time (Days)",
    censor.shape = 124,
    censor.size = 2,
    conf.int = FALSE,
    break.x.by = 720,
    ggtheme = theme_classic() +  # Use clean theme
      theme(
        axis.title.x = element_text(size = 18),  # x-axis title font size
        axis.title.y = element_text(size = 18),  # y-axis title font size
        axis.text.x = element_text(size = 18),   # x-axis tick label font size
        axis.text.y = element_text(size = 18),   # y-axis tick label font size
        legend.title = element_text(size = 18),  # Legend title
        legend.text = element_text(size = 18)
      )
  )

  # Output KM plot to PDF
  pdf(paste0(store_name, "_km.pdf"), width = 8, height = 12)
  print(lasso_km)
  dev.off()
}

# Return model and auc
svm_classfic_train <- function(seu_obj, train_score, target_class, target_prob) {
  # SVM method
  # seu_obj[[train_score]] when train_score is a column only
  train_data_x <- seu_obj[, train_score, drop = FALSE]
  # Remove NA
  valid_idx <- !is.na(seu_obj[[target_class]])
  train_data_x <- train_data_x[valid_idx, ]
  train_data_y <- factor(seu_obj[[target_class]][valid_idx])
  model <- svm(train_data_x, train_data_y, type = "C", kernel = "radial", probability = TRUE)
  pred <- predict(model, train_data_x, probability = TRUE)
  probs <- attr(pred, "probabilities")
  pred_prob <- probs[, target_prob]
  # SVM method
  roc_obj <- roc(train_data_y, pred_prob)
  train_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", train_auc))
  list(
    model = model,
    roc = roc_obj
  )
}

# Return auc
svm_classfic_predict <- function(train_model, seu_obja, train_score, test_dataset, target_class, target_prob) {
  test_auc <- 0
  valid_test_idx <- !is.na(seu_obja[[train_score]]) & !is.na(seu_obja[[target_class]]) & seu_obja[[target_class]] != ""
  seu_obja_sub <- seu_obja[valid_test_idx, ]
  if (nrow(seu_obja_sub) > 0) {
    pred <- predict(train_model, seu_obja_sub[, train_score, drop = FALSE],
      probability = TRUE
    )
    probs <- attr(pred, "probabilities")
    if (target_prob %in% colnames(probs)) {
      pred_prob <- probs[, target_prob]
      roc_obj <- roc(seu_obja_sub[[target_class]], pred_prob)
      test_auc <- as.numeric(auc(roc_obj))
      print(paste0("auc_value:",  test_auc))
    } else {
      print(paste0("Warning: 'Tumor' column not found in probabilities for ", test_dataset))
    }
  } else {
    print(paste0("Warning: No valid data for ", train_score, " in ", test_dataset))
  }
  roc_obj
}

# Return model and auc
xgboost_classfic_train <- function(seu_obj, train_score, target_class, target_prob) {
  dtrain <- xgb.DMatrix(data = as.matrix(seu_obj[, train_score, drop = FALSE]),
                        label = as.numeric(seu_obj[[target_class]] == target_prob))
  params <- list(
    objective = "binary:logistic",  # Binary classification
    eval_metric = "logloss",        # Evaluation metric
    max_depth = 6,                  # Maximum depth of trees
    eta = 0.3,                      # Learning rate
    subsample = 0.8,                # Sample sampling ratio
    colsample_bytree = 0.8          # Feature sampling ratio
  )
  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 100,                  # Number of iterations
    #early_stopping_rounds = 10,     # Early stopping
    verbose = 1
  )
  pred_prob <- predict(xgb_model, dtrain)
  roc_obj <- roc(seu_obj[[target_class]], pred_prob)
  train_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", train_auc))
  list(
    model = xgb_model,
    roc = roc_obj
  )
}

# Return auc
xgboost_classfic_predict <- function(train_model, seu_obja, train_score, test_dataset, target_class, target_prob) {
  valid_test_idx <- !is.na(seu_obja[[target_class]]) & seu_obja[[target_class]] != ""
  seu_obja_sub <- seu_obja[valid_test_idx, ]
  dtrain <- xgb.DMatrix(data = as.matrix(seu_obja_sub[, train_score, drop = FALSE], ncol = 1),
                        label = as.numeric(seu_obja_sub[[target_class]] == target_prob))
  pred_prob <- predict(train_model, dtrain)
  roc_obj <- roc(seu_obja_sub[[target_class]], pred_prob)
  test_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", test_auc))
  roc_obj
}

rf_classfic_train <- function(seu_obj, train_score, target_class, target_prob) {
  train_data <- seu_obj %>% select(all_of(c(train_score, target_class)))
  valid_idx <- !is.na(train_data[[target_class]])
  train_data <- train_data[valid_idx, ]
  train_data[[target_class]] <- as.factor(train_data[[target_class]])
  formula <- reformulate(train_score, response = target_class)
  rf_model <- randomForest(
    formula = formula,  # Dependent variable ~ All independent variables
    data = train_data,
    ntree = 500,                # Number of trees
    mtry = sqrt(ncol(train_data) - 1),  # Number of features used per tree
    importance = TRUE           # Calculate feature importance
  )
  pred_prob <- predict(rf_model, type = "prob")[, target_prob]
  roc_obj <- roc(train_data[[target_class]], pred_prob)
  train_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", train_auc))
  list(
    model = rf_model,
    roc = roc_obj
  )
}

rf_classfic_predict <- function(train_model, seu_obja, train_score, test_dataset, target_class, target_prob) {
  valid_test_idx <- !is.na(seu_obja[[target_class]]) & seu_obja[[target_class]] != ""
  seu_obja_sub <- seu_obja[valid_test_idx, ]
  test_data <- seu_obja_sub %>% select(all_of(c(train_score, target_class)))
  pred_prob <- predict(train_model, newdata = test_data, type = "prob")[, target_prob]
  roc_obj <- roc(test_data[[target_class]], pred_prob)
  test_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", test_auc))
  roc_obj
}

gbm_classfic_train <- function(seu_obj, train_score, target_class, target_prob) {
  formula <- reformulate(train_score, response = target_class)
  train_data <- seu_obj %>% select(all_of(c(train_score, target_class)))
  valid_idx <- !is.na(train_data[[target_class]])
  train_data <- train_data[valid_idx, ]
  train_data[[target_class]] <- as.numeric(train_data[[target_class]] == target_prob)
  fit <- gbm(
    formula = formula, data = train_data, distribution = "bernoulli",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 5,
    n.cores = 2
  )
  # Find index for number trees with minimum CV error
  best <- which.min(fit$cv.error)
  set.seed(123)
  gbm_model <- gbm(
    formula = formula, data = train_data, distribution = "bernoulli",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 5,
    n.cores = 2
  )
  pred_prob <- as.numeric(predict(fit, train_data, n.trees = best, type = "response"))
  roc_obj <- roc(train_data[[target_class]], pred_prob)
  train_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", train_auc))
  list(
    model = gbm_model,
    roc = roc_obj,
    best = best
  )
}

gbm_classfic_predict <- function(train_model, seu_obja, train_score, test_dataset, best, target_class, target_prob) {
  valid_test_idx <- !is.na(seu_obja[[target_class]]) & seu_obja[[target_class]] != ""
  seu_obja_sub <- seu_obja[valid_test_idx, ]
  test_data <- seu_obja_sub %>% select(all_of(c(train_score, target_class)))
  pred_prob <- as.numeric(predict(train_model, test_data, n.trees = best, type = "response"))
  roc_obj <- roc(test_data[[target_class]], pred_prob)
  test_auc <- as.numeric(auc(roc_obj))
  print(paste0("auc_value:", test_auc))
  roc_obj
}

store_path <- "store_result/"
pro_base <- "./files/"
setwd(pro_base)
# Using TCGA as SVM training model
train_scores <- list(
  c("SenMayo"),
  c("CellAge"),
  c("GenAge"),
  c("ASIG"),
  c("SASP"),
  c("AgingAtlas"),
  c("SenUp"),
  c("SigRS"),
  c("DAS"),
  c("mSS"),
  c("SenUP"),
  c("hUSI"),
  c("ClassicalSs"),
  c("HCCaging")
)
# Load other 6 samples to predict using SVM
test_datasets <- list(
  "GSE14520_sample_info_tumor_score.qs",
  "GSE76427_sample_info_tumor_score.qs",
  "OEP000321_sample_info_tumor_score.qs",
  "GSE10143_sample_info_tumor_score.qs",
  "GSE25097_sample_info_tumor_score.qs",
  "GSE63898_sample_info_tumor_score.qs",
  "LIRI_sample_info_tumor_score.qs"
)
# Empty dataframe to store the result
results_df <- data.frame()
target_class <- "sample_type"
target_prob <- "Tumor"
for (train_score in train_scores){
  print(train_score)
  cur_result <- c()
  pdf("HCCAging_ROC_RF.pdf", width = 8, height = 8)
  dataset_colors <- c(
    "TCGA" = "#f9a27d",
    "GSE14520" = "#6cb5e9",
    "GSE76427" = "#a6ce58",
    "OEP000321" = "#e48dba",
    "GSE10143" = "#bb7043",
    "GSE25097" = "#f58b4c",
    "GSE63898" = "#f9c27d",
    "LIRI" = "#a0afd3"
  )
  linetypes <- c(1, 1, 1, 1, 1, 1, 1, 1)
  model_names <- c("TCGA", "GSE14520", "GSE76427", "OEP000321", "GSE10143", "GSE25097", "GSE63898", "LIRI")

  seu_obj <- qread("taga_hcc_sample_info_all_score.qs")
  #train_result <- svm_classfic_train(seu_obj, train_score, target_class, target_prob)
  #train_result <- xgboost_classfic_train(seu_obj, train_score, target_class, target_prob)
  #train_result <- rf_classfic_train(seu_obj, train_score, target_class, target_prob)
  train_result <- gbm_classfic_train(seu_obj, train_score, target_class, target_prob)
  cur_result <- c(cur_result, as.numeric(auc(train_result$roc)))

  plot(train_result$roc,
    main = "Comparison of Four ROC Curves",
    col = dataset_colors[["TCGA"]],
    lwd = 2,
    legacy.axes = TRUE,
    print.auc = FALSE
  )
  
  for (test_dataset in test_datasets) {
    seu_obja <- qread(test_dataset)
    test_dataset_name <- toupper(strsplit(test_dataset, "_")[[1]][1])
    #test_result <- svm_classfic_predict(train_result$model, seu_obja, train_score, test_dataset, target_class, target_prob)
    #test_result <- xgboost_classfic_predict(train_result$model, seu_obja, train_score, test_dataset, target_class, target_prob)
    #test_result <- rf_classfic_predict(train_result$model, seu_obja, train_score, test_dataset, target_class, target_prob)
    test_result <- gbm_classfic_predict(train_result$model, seu_obja, train_score, test_dataset, train_result$best, target_class, target_prob)
    lines(test_result, col = dataset_colors[[test_dataset_name]], lwd = 2)
    cur_result <- c(cur_result, as.numeric(auc(test_result)))
  }
  legend("bottomright",
    legend = c(
      paste0(model_names[1], " (AUC = ", round(cur_result[1], 3), ")"),
      paste0(model_names[2], " (AUC = ", round(cur_result[2], 3), ")"),
      paste0(model_names[3], " (AUC = ", round(cur_result[3], 3), ")"),
      paste0(model_names[4], " (AUC = ", round(cur_result[4], 3), ")"),
      paste0(model_names[5], " (AUC = ", round(cur_result[5], 3), ")"),
      paste0(model_names[6], " (AUC = ", round(cur_result[6], 3), ")"),
      paste0(model_names[7], " (AUC = ", round(cur_result[7], 3), ")"),
      paste0(model_names[8], " (AUC = ", round(cur_result[8], 3), ")")
    ),
    col = dataset_colors,
    lwd = 2,
    lty = linetypes,
    bty = "n", # No border
    cex = 0.8
  )
  grid()
  dev.off()
  cur_df <- as.data.frame(t(cur_result))
  results_df <- bind_rows(results_df, cur_df)
}
rownames(results_df) <- train_scores
upper_names <- sapply(test_datasets, function(x) {
  toupper(strsplit(x, "_")[[1]][1])
})
colnames(results_df) <- c("TCGA", upper_names)

write.csv(results_df, paste0("Tumor_Normal_", "rf", ".csv"))