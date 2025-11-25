# ===== 依赖 =====
library(Seurat)
library(SeuratObject)
library(scCATCH)
library(Matrix)
library(dplyr)
library(tidyr)
library(stringr)

# 你的准确度函数（可选）
accuracy_fun <- function(pred, truth) mean(pred == truth, na.rm = TRUE)

# ---- 从 result 生成自定义 marker 表 ----
.build_custom_markers <- function(result, species = "Human", tissue = "Synthetic") {
  demo_cols <- names(scCATCH::demo_marker())
  P <- nrow(result$te); K <- ncol(result$markerGenes_error)
  gene_names <- rownames(result$te); if (is.null(gene_names)) gene_names <- paste0("Gene", seq_len(P))
  type_names <- paste0("Type", seq_len(K))
  
  zero_rows <- rowSums(result$markerGenes_error != 0) == 0
  result$markerGenes_error[zero_rows, ] <- ifelse(result$markerGenes_error[zero_rows, ] == 0, 1, 0)
  
  pairs_list <- lapply(seq_len(K), function(k) {
    idx <- which(result$markerGenes_error[, k] == 0L)  # 0 = marker
    if (length(idx) == 0) return(NULL)
    data.frame(
      species    = species,
      tissue     = tissue,
      cancer     = NA_character_,
      celltype   = type_names[k],
      subtype    = NA_character_,
      gene = gene_names[idx],
      evidence   = NA_character_,
      PMID       = NA_character_,
      source     = "Simulated",
      marker     = 1L,
      note       = NA_character_,
      stringsAsFactors = FALSE
    )
  })
  pairs <- do.call(rbind, pairs_list)
  for (miss in setdiff(demo_cols, colnames(pairs))) pairs[[miss]] <- NA
  pairs <- pairs[, demo_cols]
  pairs
}

# ---- 兜底指标 ----
.fallback_metrics <- function(pred, truth) {
  stopifnot(length(pred) == length(truth))
  acc <- mean(pred == truth, na.rm = TRUE)
  conf <- table(truth, pred, useNA = "ifany")
  classes <- sort(union(unique(truth), unique(pred)))
  prf <- lapply(classes, function(cls) {
    tp <- sum(pred == cls & truth == cls, na.rm = TRUE)
    fp <- sum(pred == cls & truth != cls, na.rm = TRUE)
    fn <- sum(pred != cls & truth == cls, na.rm = TRUE)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    recall    <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
      2 * precision * recall / (precision + recall) else NA_real_
    data.frame(class = cls, precision = precision, recall = recall, f1 = f1)
  })
  prf <- do.call(rbind, prf)
  list(accuracy = acc, confusion = conf, per_class = prf)
}

# ---- 主流程（修复：设置 VariableFeatures）----
run_scCATCH_with_accuracy <- function(
    result,
    cells_idx     = NULL,     # 可选：该组细胞列索引；不传=全部细胞
    species       = "Human",
    tissue        = "Synthetic",
    accuracy_fun  = NULL,     # 你的准确度函数 function(pred, truth) ...
    do_umap_plot  = TRUE,
    use_all_genes = TRUE,    # 若 TRUE：不用 HVG，直接把全部基因设为 VariableFeatures
    nfeatures     = 2000      # 若 use_all_genes=FALSE，则用多少 HVG
) {
  # 1) 选细胞 + 转稀疏
  idx <- if (!is.null(cells_idx)) cells_idx else seq_len(ncol(result$te))
  stopifnot(length(idx) > 0)
  expr_mat <- result$te[, idx, drop = FALSE]
  if (!inherits(expr_mat, "dgCMatrix")) expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
  
  # 2) CreateSeuratObject（不 Normalize），把矩阵放入 "data" layer（v5 正确做法）
  seu <- CreateSeuratObject(counts = expr_mat, assay = "RNA", project = "Simulated")
  DefaultAssay(seu) <- "RNA"
  seu <- SetAssayData(seu, assay = "RNA", layer = "data", new.data = expr_mat)
  
  # 3) 【关键修复】设置 VariableFeatures
  if (use_all_genes) {
    VariableFeatures(seu) <- rownames(seu)  # 直接用全部基因
  } else {
    # 在 "data" 层上选择 HVG（这不会更改你的数值，只是选子集）
    seu <- FindVariableFeatures(seu, assay = "RNA", selection.method = "vst",
                                nfeatures = nfeatures, layer = "data")
  }
  
  # 4) Scale → PCA → 邻居 → 聚类 → UMAP（默认参数）
  #    指定 features=VariableFeatures(seu) 保证不会再报 “No variable features”
  seu <- ScaleData(seu, features = VariableFeatures(seu))  # 默认读 "data"，写 "scale.data"
  seu <- RunPCA(seu, features = VariableFeatures(seu))     # 用已设定的特征
  seu <- FindNeighbors(seu)                                # 默认 dims = 1:10
  seu <- FindClusters(seu)                                 # 默认 resolution = 0.8
  # seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, graph = NULL, nn.name = NULL, features = NULL)       # 默认 dims=1:10                             # 默认 dims = 1:10
  # if (do_umap_plot) print(DimPlot(seu, reduction = "umap", label = TRUE))
  
  # 5) 自定义 marker
  custom_markers <- .build_custom_markers(result, species = species, tissue = tissue)
  
  # 6) scCATCH 注释（明确用 "data" 层）
  obj <- createscCATCH(
    data    = GetAssayData(seu, assay = "RNA", layer = "data"),
    cluster = as.character(Idents(seu))
  )
  obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = custom_markers)
  obj <- findcelltype(object = obj)
  
  # 7) 逐细胞预测与真值
  # ct_df  <- obj@celltype[, c("cluster", "cell_type")]
  # ct_map <- setNames(ct_df$cell_type, ct_df$cluster)
  # pred   <- unname(ct_map[as.character(Idents(seu))])

  
  ct <- as.data.frame(obj@celltype)
  
  col_cluster <- "cluster"
  col_type    <- "cell_type"
  col_score   <- "celltype_score"
  
  # 1) 把 “Type3, Type7, Type8” / “0.5, 0.5, 0.5” 这种列拆成长表，并把分数转数值
  ct_long <- ct %>%
    mutate(
      !!col_type := as.character(.data[[col_type]]),
      !!col_score := as.character(.data[[col_score]])
    ) %>%
    separate_rows(!!sym(col_type), !!sym(col_score), sep = ",") %>%
    mutate(
      !!col_type  := str_trim(.data[[col_type]]),
      .score_num  = as.numeric(str_trim(.data[[col_score]]))
    )
  
  # 2) 每个 cluster 取最高分；若平局，slice_max(..., with_ties = FALSE) 只保留第一个
  ct_df <- ct_long %>%
    group_by(.data[[col_cluster]]) %>%
    slice_max(order_by = .score_num, n = 1, with_ties = FALSE) %>%  # 平局→选第一个
    ungroup() %>%
    select(cluster = all_of(col_cluster),
           celltype = all_of(col_type))
  
  # 3) 映射到每个细胞（逐细胞唯一标签）
  cluster_ids <- as.character(Idents(seu))
  map_lab <- setNames(ct_df$celltype, as.character(ct_df$cluster))
  pred <- unname(map_lab[cluster_ids])
  
  # 可选：写回 meta.data
  # seu$scCATCH_single <- pred
  
  if (!is.null(result$label)) {
    truth_all <- paste0("Type", result$label)
    truth <- truth_all[idx]
  } else {
    truth <- rep(NA_character_, length(pred))
    warning("result$label 不存在，无法计算准确度；仅返回 scCATCH 结果。")
  }
  
  # 8) 准确度：优先用你传入的函数，否则兜底
  acc_result <-
    if (!is.null(accuracy_fun) && all(!is.na(truth))) {
      accuracy_fun(pred, truth)
    } else if (all(!is.na(truth))) {
      .fallback_metrics(pred, truth)
    } else {
      NA
    }
  
  # 写回 Seurat
  seu$scCATCH_pred <- pred
  
  list(
    seurat          = seu,
    scCATCH_object  = obj,
    predictions     = pred,
    truth           = truth,
    accuracy_result = acc_result
  )
}
# ====== 用法示例 ======
# 1) 用你自己的“准确度代码”（封成函数）：
# res <- run_scCATCH_with_accuracy(result, accuracy_fun = my_accuracy_fun_example)

# 2) 如果你暂时不传，就用内置兜底指标：
# res <- run_scCATCH_with_accuracy(result)

# 3) 如果你只想评估第 i 组的细胞（给出列索引）：
# idx_i <- which(some_group_vector == i)
# res_i <- run_scCATCH_with_accuracy(result, cells_idx = idx_i, accuracy_fun = my_accuracy_fun_example)

# 查看结果：
# res$accuracy_result
# res$seurat$scCATCH_pred
# DimPlot(res$seurat, group.by = "scCATCH_pred", label = TRUE)

all_metrics_catch <- list()
for (i in 1:20) {
  # 获取当前 result 数据
  result <- result_list[[i]]
  res <- run_scCATCH_with_accuracy(result)

  # 计算准确度
  true_labels <- res$truth
  annotation_labels <- res$predictions
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
  all_metrics_catch[[i]] <- list(
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


# 若未安装请先安装：install.packages(c("dplyr","purrr","ggplot2"))
library(dplyr)
library(purrr)
library(ggplot2)

# 1) 自动抓取当前环境里所有以 all_metrics_ 开头的对象（例如 all_metrics_catch、all_metrics_xx 等）
bags <- mget(ls(pattern = "^all_metrics"), inherits = TRUE)

# 2) 抽取每个方法的 accuracy，拼成长表
acc_df <- imap_dfr(bags, ~{
  # .x 是该方法的列表（长度=重复次数），.y 是对象名（方法名）
  acc <- sapply(.x, function(z) z$accuracy)   # 取出 accuracy
  data.frame(method = .y, accuracy = acc, stringsAsFactors = FALSE)
})

# 可选：美化方法名（去掉前缀 all_metrics_）
acc_df$method <- sub("^all_metrics_", "", acc_df$method)

# 3) 按中位数排序显示方法的顺序（可选）
order_by_med <- acc_df %>%
  group_by(method) %>%
  summarize(med = median(accuracy, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>% pull(method)
acc_df$method <- factor(acc_df$method, levels = order_by_med)

# 4) 画箱线图（附带散点）
p <- ggplot(acc_df, aes(x = method, y = accuracy)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.6) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Method", y = "Accuracy", title = "Accuracy across 20 runs by method") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p)

# 可选：保存图片
# ggsave("accuracy_boxplot.png", p, width = 7, height = 4, dpi = 300)