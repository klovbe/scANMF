# ===== 依赖 =====
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(scCATCH)
  library(data.table)
  library(SeuratObject)
  library(dplyr)
  library(stringr)
  library(FNN)
  library(scPred)
  library(SingleCellExperiment)
  library(scRNAseq)
  library(Seurat)
  library(ggplot2)
  library(scales)
  library(SingleR)
})

Segerstolpe <- SegerstolpePancreasData()
# 提取 label 列（细胞类型标签）
Segerstolpe_label <- colData(Segerstolpe)$'cell type'
# 获取需要去除的索引（NA 和 "unclear"）
indices_to_remove <- which(is.na(Segerstolpe_label) | Segerstolpe_label %in% c("unclassified cell", "unclassified endocrine cell", "co-expression cell"))
# 从表达矩阵中去除这些样本
Segerstolpe_filtered <- Segerstolpe[, -indices_to_remove]
# 从 Segerstolpe_label 中去除这些标签
Segerstolpe_label_filtered <- Segerstolpe_label[-indices_to_remove]

Segerstolpe_filtered <- counts(Segerstolpe_filtered)
Segerstolpe_filtered <- Segerstolpe_filtered[!duplicated(rownames(Segerstolpe_filtered)), ]
Segerstolpe_filtered <- CreateSeuratObject(Segerstolpe_filtered, project = "Segerstolpe",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
Segerstolpe_filtered[["percent.mt"]] <- PercentageFeatureSet(Segerstolpe_filtered, pattern = "^MT-")
Segerstolpe_filtered <- subset(Segerstolpe_filtered, subset = percent.mt < 20)

Segerstolpe_filtered <- NormalizeData(Segerstolpe_filtered)
Segerstolpe_filtered <- FindVariableFeatures(Segerstolpe_filtered)
Segerstolpe_filtered<- ScaleData(Segerstolpe_filtered, features = rownames(Segerstolpe_filtered))


highly_variable_genes_Segerstolpe <- VariableFeatures(Segerstolpe_filtered, n=1000)
expression_matrix_Segerstolpe <- GetAssayData(Segerstolpe_filtered, assay = "RNA")
# expression_matrix_Segerstolpe <- t(apply(expression_matrix_Segerstolpe, 1, function(x) x / (sd(x) + 1e-10)))
#scANMF——————————————————————————————————————————————————————————————————————————————————————————
checkGeneSymbols <- function(x, ...) {
  data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
}
source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
# DB file
tissue = "Pancreas" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare(db_, tissue)
gs = gs_list$gs_positive
gs$`Cancer stem cells` <- NULL
matrix_from_list <- do.call(cbind, gs)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
union_vec <- toupper(union(highly_variable_genes_Segerstolpe, unique_vector))
markers<- matrix(0, nrow = length(union_vec), ncol = length(celltypes), dimnames = list(union_vec, celltypes))
# 遍历每一列，根据基因是否为标记基因进行赋值
for (col_index in 1:ncol(matrix_from_list)) {
  for (row_index in 1:length(union_vec)) {
    gene <- union_vec[row_index]
    celltype <- celltypes[col_index]
    
    # 如果基因是标记基因，则置为1
    if (gene %in% matrix_from_list[, col_index]) {
      markers[gene, celltype] <- 1
    }
  }
}

expression_matrix_Segerstolpe<-as.matrix(expression_matrix_Segerstolpe)

# 高变基因+标记基因
#基因表达矩阵与高变基因+标记基因相交部分
common_rows <- intersect(toupper(rownames(expression_matrix_Segerstolpe)), toupper(rownames(markers)))
# common_rows 现在包含了两个矩阵中行名相同的行的标识符
common_rows1 <- rownames(expression_matrix_Segerstolpe)[toupper(rownames(expression_matrix_Segerstolpe)) %in% common_rows]
expression_matrix_common <- expression_matrix_Segerstolpe[common_rows1, , drop = FALSE]
rownames(expression_matrix_common)<-toupper(rownames(expression_matrix_common)) 
markers_common<- markers[common_rows, ]
colnames(markers_common)
#zero_cols <- which(colSums(markers_common == 0) == nrow(markers_common))
# 删除全为零的列
#markers_common<- markers_common[, -zero_cols]
mapping <- c(
  "acinar cell" = "Acinar cells",
  "alpha cell" = "Alpha cells",
  "beta cell" = "Beta cells",
  "delta cell" = "Delta cells",
  "ductal cell" = "Ductal cells",
  "epsilon cell" = "Epsilon cells",
  "gamma cell" = "Gamma (PP) cells",
  "endothelial cell" = "Endothelial cells",
  "mesenchymal" = "Mesenchymal cells",
  "PSC cell" = "Pancreatic stellate cells",
  "mast cell" = "Mast cells" 
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
Segerstolpe_label_filtered <- ifelse(Segerstolpe_label_filtered %in% names(mapping), mapping[Segerstolpe_label_filtered], Segerstolpe_label_filtered)

#——————————————————————————————————
#标签信息
unique_labels <- unique(Segerstolpe_label_filtered)
num_labels <- length(unique_labels)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells <- length(Segerstolpe_label_filtered)
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells) {
  labels <- Segerstolpe_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数

###############################################
#mutiple rounds
##################################################################################
#Laplace
# 转置：列为样本，行为特征
expression_matrix_common <- as.matrix(expression_matrix_common)
X_t <- t(expression_matrix_common)
#KNN+Gauss
k <- 300
mutual <- TRUE    # TRUE = mutual kNN, FALSE = 普通 kNN
set.seed(123)

nn <- get.knn(X_t, k = k)

gamma <- 1 / (2 * median(nn$nn.dist^2))
n <- nrow(X_t)
Q <- matrix(0, n, n)

for (i in 1:n) {
  neighbors <- nn$nn.index[i, ]
  distances <- nn$nn.dist[i, ]
  Q[i, neighbors] <- exp(-gamma * (distances^2))
}

if (mutual) {
  # Mutual kNN：仅保留互为近邻的边
  Q <- Q * t(Q)
} else {
  # 普通 kNN：对称化（并集）
  Q <- ((Q + t(Q)) / 2)
}

Q <- (Q + t(Q)) / 2   # 再次确保对称
D <- diag(rowSums(Q))
L <- D - Q             # 未归一化拉普拉斯
# 若需对称归一化，可用：
# D_inv_sqrt <- diag(1 / sqrt(diag(D) + 1e-10))
# L <- diag(n) - D_inv_sqrt %*% Q %*% D_inv_sqrt


seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
cal_scanmf <- function(alpha, beta, lambda){
result_scanmf <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Segerstolpe_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  # 应用掩码，保留的细胞标签不变，其余细胞整行变为 0
  masked_label_matrix <- label_matrix * retain_mask  # 逐行运算，FALSE 位置的整行变为 0
  colnames(masked_label_matrix)<-unique_labels
  colnames(label_matrix)<-colnames(masked_label_matrix)
  #标签和标记基因对应的细胞类型的并集
  unique_labellist<-unique(c(colnames(markers_common),colnames(label_matrix)))
  #标签以外的细胞类型
  a<- setdiff(unique_labellist, colnames(label_matrix))
  #额外的细胞类型矩阵
  label_matrix1 <- matrix(0, nrow = num_cells, ncol =length(a) )
  colnames(label_matrix1)<-a
  #总的含有标签信息细胞类型矩阵
  labelMatrix<-cbind(label_matrix1 ,masked_label_matrix )
  # 将非全为0的行中的非零元素置为0，零元素置为1
  non_zero_rows <- rowSums(labelMatrix != 0) > 0
  labelMatrix[non_zero_rows, ] <- ifelse(labelMatrix[non_zero_rows, ] == 0, 1, 0)
  #最终分析所需的标签矩阵：labelMatrix
  #标记基因对应细胞类型+标签
  b<-setdiff(unique_labellist, colnames(markers_common))
  #额外的细胞类型的标记基因矩阵
  markers_common1 <- matrix(0, nrow = nrow(markers_common), ncol =length(b) )
  colnames(markers_common1)<-b
  Markers<-cbind(markers_common1 ,markers_common)
  #找到不全为0的行
  non_zero_rows <- rowSums(Markers!= 0) > 0
  Markers[non_zero_rows, ] <- ifelse(Markers[non_zero_rows, ] == 0, 1, 0)
  #最终分析所需的标记基因矩阵：Markers
  #Markers和labelMatrix列名顺序统一
  # 按照 mat1 的列名顺序排列 mat2
  Markers<- Markers[, colnames(labelMatrix)]

  
  resultSegerstolpe<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
  cell_type_matrix_Segerstolpe  <-resultSegerstolpe[[2]]
  marker_Segerstolpe<-resultSegerstolpe[[1]]
  
  max_col_indices <- apply(cell_type_matrix_Segerstolpe, 1, which.max)
  ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
  normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix_Segerstolpe ), ncol = ncol(cell_type_matrix_Segerstolpe))
  for (j in seq_along(max_col_indices)) {
    normalized_matrix[j, max_col_indices[j]] <- 1
  }
  cell_type_matrix <- normalized_matrix
  colnames(cell_type_matrix)<-colnames(labelMatrix)
  annotation_Segerstolpe<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
  # 计算注释准确率
  true_labels<-Segerstolpe_label_filtered[!retain_mask]
  annotation_labels<-annotation_Segerstolpe[!retain_mask]
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  # 输出准确率
  print(accuracy)
  table(annotation_labels, true_labels)
  result_scanmf[[i]] <- metric_cal(true_labels, annotation_labels)
}
result_scanmf
}

#ablation
#-----------------------------------
run_ablation <- function(settings) {
  results <- list()
  for (name in names(settings)) {
    cat("Running:", name, "...\n")
    p <- settings[[name]]
    results[[name]] <- cal_scanmf(alpha = p$alpha,
                                  beta = p$beta,
                                  lambda = p$lambda)
  }
  return(results)
}



extract_metrics <- function(ablation_results) {
  df <- data.frame()
  
  for (setting in names(ablation_results)) {
    runs <- ablation_results[[setting]]
    
    for (i in seq_along(runs)) {
      res <- runs[[i]]
      
      if (is.null(res) || all(is.na(res))) {
        df <- rbind(df, data.frame(
          method = setting,
          replicate = i,
          accuracy = NA,
          weighted_f1 = NA
        ))
      } else {
        df <- rbind(df, data.frame(
          method = setting,
          replicate = i,
          accuracy = res$accuracy,
          weighted_f1 = res$weighted_f1
        ))
      }
    }
  }
  
  df$method <- factor(df$method,
                      levels=c("marker","label",
                               "marker.label","marker.laplace","label.laplace","full"))
  
  return(df)
}



ablation_settings <- list(
  marker              = list(alpha=10000, beta=0,     lambda=0),
  label               = list(alpha=0,     beta=0,     lambda=10000),
  marker.label        = list(alpha=10000, beta=0,     lambda=10000),
  marker.laplace      = list(alpha=10000, beta=10,    lambda=0),
  label.laplace       = list(alpha=0,     beta=10,    lambda=10000),
  full                = list(alpha=10000, beta=10,    lambda=10000)
)
ablation_results <- run_ablation(ablation_settings)
df_ablation <- extract_metrics(ablation_results)

name_map <- c(
  marker        = "Marker-only",
  label         = "Label-only",
  marker.label  = "Marker + Label",
  marker.laplace = "Marker + Graph",
  label.laplace = "Label + Graph",
  full          = "Full Model"
)

df_ablation$method <- name_map[df_ablation$method]

df_ablation$method <- factor(
  df_ablation$method,
  levels = c("Marker-only",
             "Label-only",
             "Marker + Label",
             "Marker + Graph",
             "Label + Graph",
             "Full Model")
)

df_long <- df_ablation %>%
  tidyr::pivot_longer(
    cols = c("accuracy", "weighted_f1"),
    names_to = "metric",
    values_to = "value"
  )

p_ablation <- ggplot(df_long, aes(x = method, y = value, color = metric)) +
  
  # 更细更专业的 boxplot
  geom_boxplot(position = position_dodge(width = 0.75),
               alpha = 0.55, width = 0.6, linewidth = 0.6,
               outlier.shape = NA) +
  
  # jitter 点更柔和
  geom_jitter(position = position_jitterdodge(
    dodge.width = 0.75,
    jitter.width = 0.12
  ),
  alpha = 0.65, size = 2) +
  
  # 专业配色（更科学、对比度佳）
  scale_color_manual(
    values = c(
      "accuracy"    = "#1F77B4",  # 蓝色（偏冷）
      "weighted_f1" = "#D62728"   # 红色（偏暖）
    ),
    labels = c("Accuracy", "Weighted F1")
  ) +
  
  # 主体主题美化（去掉网格、加粗坐标轴）
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(angle = 25, hjust = 1, size = 13),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    
    # Legend 优化
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    
    # boxplot 边距美观
    plot.margin = margin(15, 15, 15, 15)
  ) +
  
  ylab("Score")

p_ablation

p1 <- ggplot(df_ablation, aes(x=method, y=accuracy, fill=method)) +
  geom_boxplot(alpha=0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size=2, alpha=0.7) +
  scale_fill_brewer(palette="Set2") +
  theme_bw(base_size=16) +
  theme(
    axis.text.x = element_text(angle=20, hjust=1, size=13),
    legend.position = "none"
  ) +
  ylab("Accuracy") + xlab("Ablation Setting")

p1

p2 <- ggplot(df_ablation, aes(x=method, y=weighted_f1, fill=method)) +
  geom_boxplot(alpha=0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size=2, alpha=0.7) +
  scale_fill_brewer(palette="Set2") +
  theme_bw(base_size=16) +
  theme(
    axis.text.x = element_text(angle=20, hjust=1, size=13),
    legend.position = "none"
  ) +
  ylab("Weighted F1") + xlab("Ablation Setting")

p2

# 输出到 PDF（矢量图）
pdf("ablation_plot.pdf", width = 7, height = 5, useDingbats = FALSE)

print(p_ablation)

dev.off()

#####-------------------
#sensitivity
##########----------------------------------
sensitivity_alpha <- function(alpha_values,
                              beta_fixed = 10,
                              lambda_fixed = 10000) {
  results <- list()
  
  for (a in alpha_values) {
    cat("Running alpha =", a, "\n")
    res <- cal_scanmf(alpha = a,
                      beta = beta_fixed,
                      lambda = lambda_fixed)
    results[[as.character(a)]] <- res
  }
  return(results)
}
sensitivity_beta <- function(beta_values,
                             alpha_fixed = 10000,
                             lambda_fixed = 10000) {
  results <- list()
  
  for (b in beta_values) {
    cat("Running beta =", b, "\n")
    res <- cal_scanmf(alpha = alpha_fixed,
                      beta = b,
                      lambda = lambda_fixed)
    results[[as.character(b)]] <- res
  }
  return(results)
}

sensitivity_lambda <- function(lambda_values,
                               alpha_fixed = 10000,
                               beta_fixed = 10) {
  results <- list()
  
  for (l in lambda_values) {
    cat("Running lambda =", l, "\n")
    res <- cal_scanmf(alpha = alpha_fixed,
                      beta = beta_fixed,
                      lambda = l)
    results[[as.character(l)]] <- res
  }
  return(results)
}
extract_sensitivity_metrics <- function(sensitivity_list) {
  df <- data.frame()
  
  for (param in names(sensitivity_list)) {
    runs <- sensitivity_list[[param]]
    
    for (i in seq_along(runs)) {
      res <- runs[[i]]
      
      df <- rbind(df, data.frame(
        param_value = param,
        replicate = i,
        accuracy = res$accuracy,
        weighted_f1 = res$weighted_f1
      ))
    }
  }
  
  df$param_value <- as.numeric(df$param_value)
  return(df)
}

summarize_sensitivity <- function(df) {
  df %>%
    dplyr::group_by(param_value) %>%
    dplyr::summarise(
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      sd_accuracy   = sd(accuracy, na.rm = TRUE),
      mean_f1       = mean(weighted_f1, na.rm = TRUE),
      sd_f1         = sd(weighted_f1, na.rm = TRUE)
    )
}


alpha_values  <- c(100, 1000, 10000, 100000)
beta_values   <- c(0.1, 1, 10, 100)
lambda_values <- c(100, 1000, 10000, 100000)


sens_alpha <- sensitivity_alpha(alpha_values)
df_alpha <- extract_sensitivity_metrics(sens_alpha)

sens_beta  <- sensitivity_beta(beta_values)
df_beta <- extract_sensitivity_metrics(sens_beta)

sens_lambda <- sensitivity_lambda(lambda_values)
df_lambda <- extract_sensitivity_metrics(sens_lambda)

save(alpha_values, sens_alpha, beta_values, sens_beta, lambda_values, sens_lambda, file = "Ser_sens.RData")

# 提取三个 sensitivity 的 accuracy / weighted_f1
all_values <- c(df_alpha$accuracy, df_beta$accuracy, df_lambda$accuracy,
                df_alpha$weighted_f1, df_beta$weighted_f1, df_lambda$weighted_f1)

# 去掉 NA
all_values <- all_values[!is.na(all_values)]

# 自动确定上下界（不从0开始）
ymin <- min(all_values) - 0.03
ymax <- 1.0
plot_sens_line <- function(df, param_name, ylim, show_y_axis = TRUE) {
  
  df_summary <- df %>%
    group_by(param_value) %>%
    summarise(
      acc_mean = mean(accuracy, na.rm=TRUE),
      acc_sd   = sd(accuracy, na.rm=TRUE),
      f1_mean  = mean(weighted_f1, na.rm=TRUE),
      f1_sd    = sd(weighted_f1, na.rm=TRUE)
    )
  
  
  ggplot(df_summary, aes(x = param_value)) +
    # Accuracy
    geom_line(aes(y = acc_mean, color = "Accuracy"), linewidth = 0.4) +
    geom_point(aes(y = acc_mean, color = "Accuracy"), size = 1.0) +
    geom_errorbar(aes(ymin = acc_mean - acc_sd,
                      ymax = acc_mean + acc_sd,
                      color = "Accuracy"),
                  width = 0.15, linewidth = 0.3) +
    # Weighted F1
    geom_line(aes(y = f1_mean, color = "Weighted F1"), linewidth = 0.4) +
    geom_point(aes(y = f1_mean, color = "Weighted F1"), size = 1.0) +
    geom_errorbar(aes(ymin = f1_mean - f1_sd,
                      ymax = f1_mean + f1_sd,
                      color = "Weighted F1"),
                  width = 0.15, linewidth = 0.3) +
    
    scale_x_log10(breaks = df_summary$param_value) +
    scale_y_continuous(limits = ylim,
                       expand = expansion(mult = c(0.05, 0.1))) +
    
    scale_color_manual(values = c(
      "Accuracy" = "#1f78b4",
      "Weighted F1" = "#33a02c"
    )) +
    
    theme_bw(base_size = 10) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = if (show_y_axis) element_text(size=10) else element_blank(),
      axis.text.y  = if (show_y_axis) element_text(size=10) else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank()
    ) +
    xlab(param_name) +
    ylab(if (show_y_axis) "Score" else NULL)
}

p_alpha  <- plot_sens_line(df_alpha,  "Alpha",  ylim = c(ymin, ymax), show_y_axis = TRUE)
p_beta   <- plot_sens_line(df_beta,   "Gamma",   ylim = c(ymin, ymax), show_y_axis = TRUE)
p_lambda <- plot_sens_line(df_lambda, "Beta", ylim = c(ymin, ymax), show_y_axis = TRUE)

library(patchwork)
theme_set(theme_bw(base_size = 13))
theme_update(
  plot.title = element_text(size=13),
  axis.title = element_text(size=13),
  axis.text  = element_text(size=13)
)


panel_ABC <- (p_alpha | p_lambda | p_beta) +
     plot_annotation(tag_levels = "a") &
     theme(
        plot.tag = element_text(size = 10, face = "bold"),
         plot.tag.position = c(0, 1),
         plot.tag.margin = margin(t = 7)   # 把标签往下推 8 pt
      )

panel_new <- 
  ((p_alpha | p_lambda) / p_beta) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

panel_new


# 空白图用于占位
blank <- ggplot() + theme_void()

panel_final <-
  (p_alpha | p_lambda) /
  (p_beta | blank) +
  plot_layout(
    widths = c(1, 1),      # 左右宽度一致
    heights = c(1, 1)      # 上下两行高度一致（可调）
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0, 1)     # 标签左上角
  )
panel_final


# 输出到 PDF（矢量图）
pdf("sens_plot_seger.pdf", width = 7, height = 3, useDingbats = FALSE)

print(panel_ABC)

dev.off()

plot_sens_line <- function(df, param_name, ylim, show_y_axis = TRUE) {
  
  df_summary <- df %>%
    group_by(param_value) %>%
    summarise(
      acc_mean = mean(accuracy, na.rm=TRUE),
      acc_sd   = sd(accuracy, na.rm=TRUE),
      f1_mean  = mean(weighted_f1, na.rm=TRUE),
      f1_sd    = sd(weighted_f1, na.rm=TRUE)
    )
  
  ggplot(df_summary, aes(x = param_value)) +
    geom_line(aes(y = acc_mean, color = "Accuracy"), linewidth = 0.4) +
    geom_point(aes(y = acc_mean, color = "Accuracy"), size = 1.0) +
    geom_errorbar(aes(ymin = acc_mean - acc_sd,
                      ymax = acc_mean + acc_sd,
                      color = "Accuracy"),
                  width = 0.15, linewidth = 0.3) +
    
    geom_line(aes(y = f1_mean, color = "Weighted F1"), linewidth = 0.4) +
    geom_point(aes(y = f1_mean, color = "Weighted F1"), size = 1.0) +
    geom_errorbar(aes(ymin = f1_mean - f1_sd,
                      ymax = f1_mean + f1_sd,
                      color = "Weighted F1"),
                  width = 0.15, linewidth = 0.3) +
    
    scale_x_log10(breaks = df_summary$param_value) +
    scale_y_continuous(limits = ylim,
                       expand = expansion(mult = c(0.05, 0.1))) +
    
    scale_color_manual(values = c(
      "Accuracy" = "#1f78b4",
      "Weighted F1" = "#33a02c"
    )) +
    
    theme_bw(base_size = 13) +
    theme(
      legend.position = c(0.02, 0.02),   # 图内左下角
      legend.justification = c(0, 0),
      legend.title = element_blank(),
      legend.background = element_rect(fill = alpha("white", 0.5), color = NA),
      legend.key.size = unit(0.6, "lines"),
      
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.y = if (show_y_axis) element_text(size=13) else element_blank(),
      axis.text.y  = if (show_y_axis) element_text(size=13) else element_blank(),
      axis.ticks.y = if (show_y_axis) element_line() else element_blank()
    ) +
    xlab(param_name) +
    ylab(if (show_y_axis) "Score" else NULL)
}

p_alpha  <- plot_sens_line(df_alpha,  "Alpha",  ylim = c(ymin, ymax), show_y_axis = TRUE)
p_beta   <- plot_sens_line(df_beta,   "Gamma",   ylim = c(ymin, ymax), show_y_axis = TRUE)
p_lambda <- plot_sens_line(df_lambda, "Beta", ylim = c(ymin, ymax), show_y_axis = TRUE)



ggsave(
  filename = "sensalpha.png",   # 你的文件名
  plot     = p_alpha,       # 你的 ggplot 对象
  width    = 7,                 # 图的宽度（单位：英寸）
  height   = 6,                 # 图的高度（单位：英寸）
  units    = "in",              # 或 "cm"
  dpi      = 300,               # 出版级分辨率（300–600）
  bg       = "white"            # 背景色（避免透明导致变灰）
)
ggsave(
  filename = "sensgamma.png",   # 你的文件名
  plot     = p_beta,       # 你的 ggplot 对象
  width    = 7,                 # 图的宽度（单位：英寸）
  height   = 6,                 # 图的高度（单位：英寸）
  units    = "in",              # 或 "cm"
  dpi      = 300,               # 出版级分辨率（300–600）
  bg       = "white"            # 背景色（避免透明导致变灰）
)

ggsave(
  filename = "sensbeta.png",   # 你的文件名
  plot     = p_lambda,       # 你的 ggplot 对象
  width    = 7,                 # 图的宽度（单位：英寸）
  height   = 6,                 # 图的高度（单位：英寸）
  units    = "in",              # 或 "cm"
  dpi      = 300,               # 出版级分辨率（300–600）
  bg       = "white"            # 背景色（避免透明导致变灰）
)


