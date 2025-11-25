#模拟数据生成算法
#K-细胞类型数，Kn_vec-每类细胞个数向量，Ndiff-差异表达基因数，Nsame-相似表达基因数,logMean, 正态分布参数，type="cluster", "DE"or "marker"，
# ZeroRate-0的比例，a-保留标签比例，b-已知标签中引入错误的比例，errpr_rate-标记错误比例，sigmahetero控制噪声随表达强度变化的异方差部分，sigmahomo是恒定噪声项

model_sample <- function(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero,sigmahomo,type, a, b, error_rate) {
  K <- K
  N <- sum(Kn_vec)  # Total number of cells
  P <- Ndiff + Nsame  # Total number of genes
  
  # Generate the labels with an uneven distribution based on Kn_vec
  Label <- unlist(mapply(function(k, n) rep(k, n), 1:K, Kn_vec))  # Uneven labels
  
  Z <- matrix(0, N, K)  # Initialize Z matrix
  start <- 1
  for (i in 1:K) {
    end <- start + Kn_vec[i] - 1
    Z[start:end, i] <- 1  # Assign Z values for each cell type
    start <- end + 1
  }
  
  # Generate expression values for genes
  Esame <- rnorm(Nsame, logMean, logSd)
  Esame <- Esame * rbinom(Nsame, 1, 1 - ZeroRate)  # Simulate true zeros
  Esame <- rep(Esame, K)
  dim(Esame) <- c(Nsame, K)
  
  if (type == "cluster") {
    Ediff <- matrix(rnorm(K * Ndiff, logMean, logSd), nrow = K, ncol = Ndiff)
  } else if (type == "DE") {
    Ediff <- rnorm(Ndiff, logMean, logSd)
    Ediff <- matrix(Ediff, nrow = K, ncol = Ndiff, byrow = TRUE)
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] <- Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] * 2
    }
  } else if (type == "marker") {
    Ediff <- matrix(rnorm(K * Ndiff, logMean - 0.5, logSd), nrow = K, ncol = Ndiff)
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] <- rnorm(Ndiff.u, logMean, logSd) * 1.5
    }
  }
  Ediff <- Ediff * matrix(rbinom(K * Ndiff, 1, 1 - ZeroRate), K, Ndiff)
  
  Te <- Z %*% cbind(Ediff, t(Esame))
  Te[Te > 0] <- Te[Te > 0] + rnorm(sum(Te > 0), 0, Te[Te > 0] * sigmahetero + sigmahomo)
  Te[Te <= 0] <- rnorm(sum(Te <= 0), 0, 1) - 1.5
  Te[Te <= 0] <- 0
  
  De <- t(Te)
  
  # Marker gene matrix
  markerGenes <- matrix(1, nrow = K, ncol = P)  # Initialize as 1
  if (type == "marker") {
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      for (j in 1:Ndiff.u) {
        markerGenes[k, ((k - 1) * Ndiff.u + j)] <- 0  # Set marker gene positions to 0
      }
    }
  }
  
  # Generate sample, type, and gene names
  sample_names <- paste("Sample", 1:N, sep = "")
  type_names <- paste("Type", 1:K, sep = "")
  gene_names <- paste("Gene", 1:P, sep = "")
  
  rownames(markerGenes) <- type_names
  colnames(markerGenes) <- gene_names
  markerGenes <- t(markerGenes)
  rows_all_ones <- apply(markerGenes, 1, function(x) all(x == 1))
  markerGenes[rows_all_ones, ] <- 0
  
  rownames(Te) <- sample_names
  colnames(Te) <- gene_names
  
  # Marker gene information
  marker_gene_info <- list()
  for (i in 1:P) {
    marker_in_types <- which(markerGenes[i, ] == 0)
    if (length(marker_in_types) > 0) {
      marker_gene_info[[gene_names[i]]] <- type_names[marker_in_types]
    }
  }
  
  # Transpose gene expression matrix
  Te <- t(Te)
  
  # Generate correct label matrix
  num_cell_types <- length(unique(Label))
  label_matrix_correct <- matrix(0, nrow = N, ncol = num_cell_types)
  for (i in 1:num_cell_types) {
    label_matrix_correct[Label == i, i] <- 1
  }
  
  # Generate label matrix with errors
  label_matrix_with_errors <- label_matrix_correct
  percentage <- a
  for (i in 1:num_cell_types) {
    cell_indices <- which(Label == i)
    num_to_keep <- ceiling(percentage * length(cell_indices))
    keep_indices <- sample(cell_indices, num_to_keep)
    drop_indices <- setdiff(cell_indices, keep_indices)
    label_matrix_with_errors[drop_indices, ] <- 0
    label_matrix_with_errors[keep_indices, i] <- 1
  }
  
  for (i in 1:num_cell_types) {
    cell_indices <- which(label_matrix_with_errors[, i] == 1)
    num_errors <- ceiling(b * length(cell_indices))
    error_indices <- sample(cell_indices, num_errors)
    for (idx in error_indices) {
      other_types <- setdiff(1:num_cell_types, i)
      label_matrix_with_errors[idx, i] <- 0
      label_matrix_with_errors[idx, sample(other_types, 1)] <- 1
    }
  }
  
  celltype_matrix <- label_matrix_with_errors
  for (i in 1:N) {
    if (all(celltype_matrix[i, ] == 0)) {
      celltype_matrix[i, ] <- label_matrix_correct[i, ]
    }
  }
  
  # Generate marker genes with errors
  generate_marker_with_errors <- function(n_genes, n_cell_types, error_rate) {
    total_elements <- n_genes * n_cell_types
    num_errors <- round(total_elements * error_rate)
    
    error_indices <- sample(1:total_elements, num_errors)
    markerGenes_error <- markerGenes
    for (index in error_indices) {
      row_index <- ceiling(index / n_cell_types)
      col_index <- index %% n_cell_types
      if (col_index == 0) col_index <- n_cell_types
      markerGenes_error[row_index, col_index] <- 1 - markerGenes_error[row_index, col_index]
    }
    return(markerGenes_error)
  }
  
  markerGenes_error <- generate_marker_with_errors(P, K, error_rate)
  
  # Return simulation data
  sData <- list(K = K, N = N, P = P, label = Label, te = Te, de = De, markerGenes_error = markerGenes_error,
                labelMatrixCorrect = label_matrix_correct, labelMatrixWithErrors = label_matrix_with_errors,
                celltype_matrix = celltype_matrix, marker_gene_info = marker_gene_info)
  return(sData)
}
# -----------------------------------------------------------
# 初始化一个空列表来存储模拟数据生成结果
result_list <- list()
# 设置参数
K <- 8# 细胞类型数
Kn_vec<-c(300,600,900,200,100,300,100,50)#每个类别细胞数
Ndiff <- 32 # 差异表达基因数
Nsame <- 400 # 相同表达基因数
logMean <- 1 # 均值
logSd <- 0.1 # 标准差
ZeroRate <- 0.1 # 零值比例
type <- "marker" # 差异表达类型
sigmahetero<-0.2
sigmahomo<-0.1
a <- 0.1 # 标签保留比例
# 生成20组数据

#一个批次
for (i in 1:20) {
  # model_sample <- function(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero,sigmahomo,type, rate to keep the labels, label error, marker error_rate) 
  result <- model_sample(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero,sigmahomo,type, a, 0, 0)
  result_list[[i]] <- result
}

# 查看生成的结果
print(length(result_list))  # 输出生成的结果组数，应该是20
print(str(result_list[[1]]))  # 查看第1组数据的结构
#result_moni1<-result_list[[1]]
#result_moni2<-result_list[[2]]
#result_moni3<-result_list[[3]]
#result_moni4<-result_list[[4]]
#result_moni5<-result_list[[5]]
#result_moni6<-result_list[[6]]
#result_moni7<-result_list[[7]]
#result_moni8<-result_list[[8]]
#result_moni9<-result_list[[9]]
#result_moni10<-result_list[[10]]
#result_moni11<-result_list[[11]]
#result_moni12<-result_list[[12]]
#result_moni13<-result_list[[13]]
#result_moni14<-result_list[[14]]
#result_moni15<-result_list[[15]]
#result_moni16<-result_list[[16]]
#result_moni17<-result_list[[17]]
#result_moni18<-result_list[[18]]
#result_moni19<-result_list[[19]]
#result_moni20<-result_list[[20]]

# 创建一个空列表来保存每组的结果
all_metrics_scANMF <- list()

# 遍历 result_list 中的每个 result
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, t(labelMatrix), ncol(labelMatrix), markerGenes, 1000, 0.04, 1200, L, D, Q, 100, 0.1)
  # result_moni <- scANMF(expression, t(labelMatrix), ncol(labelMatrix), markerGenes, 10, 0.04, 1200, L, D, Q, 1000, 0.1)
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics_scANMF[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}

------------------------------------------------------------
#参数敏感性分析，变换参数值组合，画图

# 提取每一组的 accuracy 值
accuracies0.05 <- sapply(all_metrics, function(x) x$accuracy)
weightedf10.05 <- sapply(all_metrics, function(x) x$weighted_f1)


accuracies800 <- sapply(all_metrics, function(x) x$accuracy)
weightedf1800 <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies600 <- sapply(all_metrics, function(x) x$accuracy)
weightedf1600 <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1200 <- sapply(all_metrics, function(x) x$accuracy)
weightedf11200 <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1400 <- sapply(all_metrics, function(x) x$accuracy)
weightedf11400 <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1800l <- sapply(all_metrics, function(x) x$accuracy)
weightedf11800l <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1400l <- sapply(all_metrics, function(x) x$accuracy)
weightedf11400l <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1200l <- sapply(all_metrics, function(x) x$accuracy)
weightedf11200l <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies1000l <- sapply(all_metrics, function(x) x$accuracy)
weightedf11000l <- sapply(all_metrics, function(x) x$weighted_f1)


accuracies1600l <- sapply(all_metrics, function(x) x$accuracy)
weightedf11600l <- sapply(all_metrics, function(x) x$weighted_f1)

accuracies<-cbind(accuracies0.02,accuracies0.03,accuracies0.04,accuracies0.05,accuracies0.06)
accuracies1<-cbind(accuracies600,accuracies800,accuracies0.04,accuracies1200,accuracies1400)
accuracies2<-cbind(accuracies1800l,accuracies1000l,accuracies1200l,accuracies1400l,accuracies1600l )



accuracy_means2 <- colMeans(accuracies2)
accuracy_sds2 <- apply(accuracies2, 2, sd)

accuracy_means <- colMeans(accuracies)
accuracy_sds <- apply(accuracies, 2, sd)

accuracy_means1 <- colMeans(accuracies1)
accuracy_sds1 <- apply(accuracies1, 2, sd)

accuracy_means2 <- colMeans(accuracies2)
accuracy_sds2 <- apply(accuracies2, 2, sd)


# Step 2: 准备横轴标签（列名）
x_labels <- c(0.02,0.03,0.04,0.05,0.06)
x_labels1 <- c(600,800,1000,1200,1400)
x_labels2 <- c(800,1000,1200,1400,1600)

# Step 3: 绘图
library(ggplot2)

# 整理成 data.frame 方便 ggplot2 使用
df <- data.frame(
  Condition = factor(x_labels2, levels = x_labels2),
  Mean = accuracy_means,
  SD = accuracy_sds
)
plot(
  x = x_labels2,
  y = accuracy_means2,
  type = "o",
  col = "darkblue",
  pch = 16,
  lwd = 2,
  xlab = expression(lambda),
  ylab = "Mean Accuracy",
  ylim = c(0.9, 1.00)
)
arrows(
  x0 = x_labels2,
  y0 = accuracy_means2 - accuracy_sds2,
  x1 = x_labels2,
  y1 = accuracy_means2 + accuracy_sds2,
  angle = 90,
  code = 3,
  length = 0.05,
  col = "gray40"
)
# Step 4: 绘制带误差条的折线图
# 继续使用前面创建的 df 数据框
library(ggplot2)

ggplot(df, aes(x = Condition, y = Mean, group = 1)) +
  geom_line(color = "#1f77b4", size = 1.3) +
  geom_point(size = 3, shape = 21, fill = "white", color = "#1f77b4", stroke = 1.2) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.15, color = "#999999", linewidth = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  ) +
  labs(
    x = expression(alpha),
    y = "Mean Accuracy"
  ) +
  ylim(0.7, 1.00)  # 
# 查看所有 accuracy
print(accuracies)
weighted_f1scores<-sapply(all_metrics, function(x) x$weighted_f1)
print(weighted_f1scores)


--------------------------------------
#跨批次

result_batch_list <- list()  # 存储结果

# 设置批次效应的强度
batch_effect_strength <- 0.2

for (i in 1:20) { 
  # 获取原始数据 
  result <- result_list[[i]] 
  expression <- result[["te"]] 
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]] 
  
  # 找到不全为0的行索引给表达值添加批次效应 
  non_zero_rows <- which(rowSums(labelMatrixWithErrors != 0) > 0)  # 获取不全为零的行的索引
  
  # 创建批次标签
  N <- nrow(labelMatrixWithErrors)  # 获取样本数量
  batch_labels <- rep(2, N)  # 默认所有样本为批次2
  batch_labels[non_zero_rows] <- 1  # 不全为零的行分配为批次1
  
  # 为不全为零的样本分配批次效应（例如，批次1和批次2）
  batch_offsets <- rnorm(2, mean = 0, sd = batch_effect_strength)  # 为两个批次生成偏移量
  
  # 将批次效应加入表达数据
  Te_with_batch_effect <- result$te
  for (row_idx in non_zero_rows) {  # 遍历不全为零的行索引
    Te_with_batch_effect[,row_idx] <- Te_with_batch_effect[, row_idx] + batch_offsets[batch_labels[row_idx]]
  }
  
  # 将修改后的数据存入结果
  result$te_with_batch_effect <- Te_with_batch_effect
  result$batch <- batch_labels  # 将批次标签存储在结果中
  
  # 将每组结果存储到 result_batch_list 中
  result_batch_list[[i]] <- result
}

result_batch1<-result_batch_list[[1]]
result_batch2<-result_batch_list[[2]]
result_batch3<-result_batch_list[[3]]
result_batch4<-result_batch_list[[4]]
result_batch5<-result_batch_list[[5]]
result_batch6<-result_batch_list[[6]]
result_batch7<-result_batch_list[[7]]
result_batch8<-result_batch_list[[8]]
result_batch9<-result_batch_list[[9]]
result_batch10<-result_batch_list[[10]]
result_batch11<-result_batch_list[[11]]
result_batch12<-result_batch_list[[12]]
result_batch13<-result_batch_list[[13]]
result_batch14<-result_batch_list[[14]]
result_batch15<-result_batch_list[[15]]
result_batch16<-result_batch_list[[16]]
result_batch17<-result_batch_list[[17]]
result_batch18<-result_batch_list[[18]]
result_batch19<-result_batch_list[[19]]
result_batch20<-result_batch_list[[20]]



all_metrics1<-list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_batch_list[[i]]
  
  # 提取相关数据
  expression <- result[["te_with_batch_effect"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 100, 0.01, 100, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics1[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}
#结果分析同一个批次

----------------------------------------
#类别均衡---消融实验、

# 初始化一个空列表
result_list2 <- list()

# 设置参数
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
a <- 0.1 # 标签保留比例
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



b_values <- c(0, 0.2, 0.4, 0.6)
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


#消融实验
#######################

result_list555 <- list()
for (i in 1:20) {
  result <- model_sample(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero, sigmahomo, type, a, 0, 0)
  result_list555[[i]] <- result
}
all_metrics1 <- list()
# 遍历 result_list 中的每个 result
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 0, 1, 100, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics1[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#没有alpha
all_metrics2 <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 100, 0, 10, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics2[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#没有beta
all_metrics3 <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 100, 1, 0, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics3[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#没有lambda
all_metrics4 <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 100, 0, 0, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics4[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#没有lambda、beta
all_metrics5 <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 0, 1, 100, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics5[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#没有alpha、beta
all_metrics6 <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list555[[i]]
  # 提取相关数据
  expression <- result[["te"]]
  markerGenes <- result[["markerGenes_error"]]
  labelMatrixWithErrors <- result[["labelMatrixWithErrors"]]
  
  # 处理 labelMatrixWithErrors，非零值行改为 1
  non_zero_rows <- rowSums(labelMatrixWithErrors != 0) > 0
  labelMatrixWithErrors[non_zero_rows, ] <- ifelse(labelMatrixWithErrors[non_zero_rows, ] == 0, 1, 0)
  labelMatrix <- labelMatrixWithErrors
  
  # 计算欧几里得距离
  distance_matrix <- as.matrix(dist(t(expression), method = "euclidean"))
  
  # 计算 gamma
  gamma <- 1 / (2 * median(distance_matrix^2))
  
  # 计算 Q 矩阵
  Q <- exp(-gamma * distance_matrix^2)
  
  # 计算度矩阵 D 和拉普拉斯矩阵 L
  D <- diag(rowSums(Q))
  L <- D - Q
  
  # 处理 label 信息
  labelmoni <- paste0("Type", result[["label"]])
  
  # 调用 YNMF 函数
  result_moni <- scANMF(expression, labelMatrix, ncol(labelMatrix), markerGenes, 1000, 1, 100, L, D, Q, 500, 0.1)
  
  # 提取结果
  cell_type_matrix_moni <- result_moni[[2]]
  colnames(cell_type_matrix_moni) <- colnames(markerGenes)
  # 计算预测标签
  max_col_names <- colnames(cell_type_matrix_moni)[apply(cell_type_matrix_moni, 1, which.max)]
  
  # 计算准确度
  true_labels <- labelmoni
  annotation_labels <- max_col_names
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  # 获取所有可能的类别标签
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  
  # 确保标签都是 factor 且 levels 对齐
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
  # 混淆矩阵
  conf_matrix <- table(Predicted = annotation_labels, Actual = true_labels)
  
  # 提取 TP, FP, FN
  TP <- diag(conf_matrix)
  FP <- rowSums(conf_matrix) - TP
  FN <- colSums(conf_matrix) - TP
  
  # 精确率和召回率
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # F1 分数
  f1_score <- 2 * precision * recall / (precision + recall)
  f1_score[is.na(f1_score)] <- 0
  
  # micro average
  mean_f1_score <- mean(f1_score)
  
  # weighted F1
  support <- colSums(conf_matrix)  # 每类真实样本数量
  total <- sum(support)
  weighted_f1 <- sum(f1_score * support) / total
  
  # 存储每组结果的指标
  all_metrics6[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  # 打印每组的结果
  cat("Results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}#都有
