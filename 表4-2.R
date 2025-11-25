K <- 4# 细胞类型数
Kn_vec<-c(300,300,300,300)
Ndiff <- 20 # 差异表达基因数
Nsame <- 100 # 相同表达基因数
logMean <- 1 # 均值
logSd <- 0.5 # 标准差
ZeroRate <- 0.1 # 零值比例
type <- "marker" # 差异表达类型
sigmahetero<-0.2
sigmahomo<-0.1
a <- 0.1 # 标签保留比例。后续换成0.2，0.3，重复下列步骤
b_values<-c(0,0.2,0.4,0.6)
error_rate_values<-c(0.2,0.4,0.6)
# 生成10组数据
result_list2_all <- list()
for (b in b_values) {
  result_list2 <- list()
  for (i in 1:10) {
    result <- model_sample(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero, sigmahomo, type, a, b, 0)
    result_list2[[i]] <- result
  }
  result_list2_all[[as.character(b)]] <- result_list2
}
for (error_rate in error_rate_values) {
  result_list3 <- list()
  for (i in 1:10) {
    result <- model_sample(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero, sigmahomo, type, a, 0, error_rate)
    result_list2[[i]] <- result
  }
  result_list3_all[[as.character(error_rate)]] <- result_list3
}
all_metrics_all <- list()  # 用于保存所有b的评估指标
for (b in b_values) {
  all_metrics <- list()  # 当前b值下的10组指标
  
  for (i in 1:10) {
    result <- result_list2_all[[as.character(b)]][[i]]
    
    expression <- result[["te"]]
    markerGenes <- result[["markerGenes_error"]]
    labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
    
    # 修正 labelMatrix
    non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
    labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
    labelMatrix <- labelMatrixWithErrors
    
    # 欧几里得距离和Laplacian计算
    distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
    gamma <- 1 / (2 * median(distance_matrix^2))
    Q <- exp(-gamma * distance_matrix^2)
    D <- diag(rowSums(Q))
    L <- D - Q
    
    # 真标签字符串
    labelmoni <- paste0("Type", result[["label"]])
    
    # 调用模型
    result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 10, 0.1, 100, L, D, Q, 500, 0.1)
    cell_type_matrix_moni <- result_moni[[2]]
    colnames(cell_type_matrix_moni) <- colnames(markerGenes)
    
    # 预测标签
    max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
    
    # 准确度计算
    true_labels <- factor(labelmoni)
    annotation_labels <- factor(max_col_names, levels = levels(true_labels))
    
    accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
    
    # 混淆矩阵
    conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
    TP <- diag(conf_matrix)
    FP <- rowSums(conf_matrix) - TP
    FN <- colSums(conf_matrix) - TP
    
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    f1_score <- 2 * precision * recall / (precision + recall)
    f1_score[is.na(f1_score)] <- 0
    
    mean_f1_score <- mean(f1_score)
    support <- colSums(conf_matrix)
    total <- sum(support)
    weighted_f1 <- sum(f1_score * support) / total
    
    # 存储每组结果
    all_metrics[[i]] <- list(
      accuracy = accuracy,
      macro_f1 = mean_f1_score,
      weighted_f1 = weighted_f1,
      f1_score_per_class = f1_score
    )
    
    cat("Results for b =", b, ", group", i, ":\n")
    cat("Accuracy:", accuracy, "\n")
    cat("Macro F1:", mean_f1_score, "\n")
    cat("Weighted F1:", weighted_f1, "\n\n")
  }
  
  # 存入总结果
  all_metrics_all[[as.character(b)]] <- all_metrics
}
#平均指标
for (b in names(all_metrics_all)) {
  accuracies <- sapply(all_metrics_all[[b]], function(x) x$accuracy)
  macro_f1s <- sapply(all_metrics_all[[b]], function(x) x$macro_f1)
  weighted_f1s <- sapply(all_metrics_all[[b]], function(x) x$weighted_f1)
  
  cat("Summary for b =", b, "\n")
  cat("Mean Accuracy:", mean(accuracies), "\n")
  cat("Mean Macro F1:", mean(macro_f1s), "\n")
  cat("Mean Weighted F1:", mean(weighted_f1s), "\n\n")
}
#error_rate相同
