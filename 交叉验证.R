# run_model.R
#数据准备
expression_matrix<-read.csv("GSE122108/expression_matrix_common.csv",row.names = 1,header = TRUE)
labelMatrix<-read.csv("GSE122108/labelMatrix.csv",row.names = 1,header=TRUE)
Markers<-read.csv("GSE122108/Markers.csv",row.names = 1,header =TRUE)
Labels_Matrix<-read.csv("GSE122108/Labels_Matrix.csv",row.names = 1,header =TRUE)
expression_matrix<-as.matrix(expression_matrix)
labelMatrix<-as.matrix(labelMatrix)
Markers<-as.matrix(Markers)
Labels_Matrix<-as.matrix(Labels_Matrix)

# 定义欧几里得距离函数
euclidean_dist <- function(matrix) {
  dist_matrix <- as.matrix(dist(matrix, method = "euclidean"))
  return(dist_matrix)
}
#粗粒度搜索空间，后续根据结果细化修改搜索空间
alpha_values <- c(0.01, 0.1, 1, 10, 100, 1000)
beta_values <- c(0.01, 0.1, 1, 10, 100, 1000)
lambda_values <- c(0.01, 0.1, 1, 10, 100, 1000)

# 初始化存储最佳结果的变量
best_alpha <- NULL
best_beta <- NULL
best_lambda <- NULL
best_coefficient <- 0

# 初始化结果数据框
results <- data.frame(
  alpha = numeric(length(alpha_values) * length(beta_values) * length(lambda_values)),
  beta = numeric(length(alpha_values) * length(beta_values) * length(lambda_values)),
  lambda = numeric(length(alpha_values) * length(beta_values) * length(lambda_values)),
  coefficient = numeric(length(alpha_values) * length(beta_values) * length(lambda_values))
)

counter <- 1

# 筛选出有已知标签信息的样本
labeled_indices <- which(rowSums(labelMatrix) > 0)
labeled_expression_matrix <- expression_matrix[, labeled_indices]
labeled_labelMatrix <- labelMatrix[labeled_indices, ]
labeled_LabelsMatrix <- Labels_Matrix[labeled_indices, ]

# 设置交叉验证
set.seed(123)
k <- 5
shuffled_indices <- sample(seq_along(labeled_indices))
fold_size <- floor(length(labeled_indices) / k)
folds <- vector("list", k)
for (i in 1:k) {
  if (i < k) {
    folds[[i]] <- shuffled_indices[((i - 1) * fold_size + 1):(i * fold_size)]
  } else {
    folds[[i]] <- shuffled_indices[((i - 1) * fold_size + 1):length(shuffled_indices)]
  }
}

# 构造拉普拉斯矩阵
s <- 1
labeled_Q <- exp(-euclidean_dist(t(labeled_expression_matrix))^2 / s^2)
labeled_D <- diag(rowSums(labeled_Q))
labeled_L <- labeled_D - labeled_Q

# 网格搜索
for (alpha in alpha_values) {
  for (beta in beta_values) {
    for (lambda in lambda_values) {
      
      coefficients <- c()
      coefficients1 <- c()
      
      for (fold in folds) {
        train_indices <- setdiff(seq_along(labeled_indices), fold)
        test_indices <- fold
        
        test_labelMatrix <- labeled_labelMatrix
        test_labelMatrix[test_indices, ] <- 0  # 屏蔽测试集标签
        
        # 执行模型训练
        test_result <- MLLNMF(
          labeled_expression_matrix, test_labelMatrix, 
          ncol(test_labelMatrix), Markers, 
          alpha, beta, lambda, 
          labeled_L, labeled_D, labeled_Q, 
          maxIter = 6000, tol = 0.00001
        )
        
        cell_type_matrix <- test_result[[2]]
        
        # 取最大值赋1，其他为0（one-hot）
        max_col_indices_test <- apply(cell_type_matrix, 1, which.max)
        normalized_matrix_test <- matrix(0, nrow = nrow(cell_type_matrix), ncol = ncol(cell_type_matrix))
        for (j in seq_along(max_col_indices_test)) {
          normalized_matrix_test[j, max_col_indices_test[j]] <- 1
        }
        
        cell_type_matrix <- normalized_matrix_test
        test_cell_type_matrix <- cell_type_matrix[test_indices, ]
        test_LabelsMatrix <- labeled_LabelsMatrix[test_indices, ]
        
        # 全部样本准确率
        comparison <- +(apply(cell_type_matrix == labeled_LabelsMatrix, 1, all))
        coefficient <- sum(comparison) / nrow(cell_type_matrix)
        coefficients <- c(coefficients, coefficient)
        
        # 验证集样本准确率
        comparison1 <- +(apply(test_cell_type_matrix == test_LabelsMatrix, 1, all))
        coefficient1 <- sum(comparison1) / nrow(test_cell_type_matrix)
        coefficients1 <- c(coefficients1, coefficient1)
      }
      
      mean_coefficient <- mean(coefficients, na.rm = TRUE)
      mean_coefficient1 <- mean(coefficients1, na.rm = TRUE)
      
      results[counter, ] <- c(alpha, beta, lambda, mean_coefficient)
      counter <- counter + 1
      
      if (mean_coefficient1 > best_coefficient) {
        best_alpha <- alpha
        best_beta <- beta
        best_lambda <- lambda
        best_coefficient <- mean_coefficient1
      }
    }
  }
}

# 输出最佳参数组合及其对应的coefficient值
cat("最佳参数组合:\n")
cat("Alpha:", alpha, "\n")
cat("Beta:", best_beta, "\n")
cat("Lambda:", best_lambda, "\n")
cat("Coefficient:", best_coefficient, "\n")

# 输出所有参数组合及其对应的coefficient值
output_file <- paste0("results_alpha_", alpha, ".csv")
write.csv(results, output_file, row.names = FALSE)
