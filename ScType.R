# =========================
# scType：与 scCATCH 同风格主流程
# 0/1 矩阵：0 = marker
# 输出字段与 run_scCATCH_with_accuracy 对齐：
#   - seurat
#   - predictions            （逐细胞唯一标签）
#   - predictions_cluster    （按 cluster 汇总后扩散到细胞）
#   - truth
#   - accuracy_result
#   - 其余：score_matrix / cluster_table / db_df / db_path
# =========================
checkGeneSymbols <- function(x, ...) {
  data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
}

run_scType_with_accuracy <- function(
    result,
    cells_idx     = NULL,     # 可选：该组细胞列索引；不传=全部细胞
    tissue        = "Custom", # 自定义库用到的 tissue 名
    out_name_cell = "ScType_single",
    out_name_clu  = "ScType_cluster",
    min_score     = NULL,     # 可选：逐细胞最大分的阈值，小于阈值标 "Unknown"
    use_all_genes = TRUE,     # TRUE：全部基因设为 VariableFeatures；与 scCATCH 参数一致
    nfeatures     = 2000,     # use_all_genes=FALSE 时使用
    cluster_method = "sum",  # 聚合方式 c("sum","mean","median")
    accuracy_fun  = NULL      # 可选：function(pred, truth) -> 指标/数值
){
  # ---- 依赖 ----
  pkgs <- c("Seurat","SeuratObject","Matrix","dplyr","openxlsx","stringr")
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install)) install.packages(to_install)
  lapply(pkgs, require, character.only = TRUE)
  
  cluster_method <- match.arg(cluster_method)
  
  # ---- 输入检查与抽子集 ----
  if (is.null(result$te)) stop("result$te 缺失：需要基因×细胞的表达矩阵")
  if (is.null(result$markerGenes_error)) stop("result$markerGenes_error 缺失：需要0/1矩阵（0=marker）")
  if (is.null(rownames(result$markerGenes_error))) stop("markerGenes_error 需要行名=基因")
  if (is.null(colnames(result$markerGenes_error))) stop("markerGenes_error 需要列名=细胞类型")
  
  idx <- if (!is.null(cells_idx)) cells_idx else seq_len(ncol(result$te))
  stopifnot(length(idx) > 0)
  
  expr_mat <- result$te[, idx, drop = FALSE]
  if (!inherits(expr_mat, "dgCMatrix")) expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
  
  # ---- 构建 scType 自定义库（0=正向 marker；负向留空）----
  zero_rows <- rowSums(result$markerGenes_error != 0) == 0
  result$markerGenes_error[zero_rows, ] <- ifelse(result$markerGenes_error[zero_rows, ] == 0, 1, 0)
  M <- as.data.frame(result$markerGenes_error, check.names = FALSE)
  for (j in seq_len(ncol(M))) M[[j]] <- as.numeric(as.character(M[[j]]))
  gene_names <- toupper(rownames(M))
  cell_types <- colnames(M)
  
  db_df <- dplyr::bind_rows(lapply(cell_types, function(ct){
    pos_genes <- gene_names[ M[[ct]] == 0 ]  # 0=marker
    dplyr::tibble(
      tissueType = tissue,
      cellName = ct,
      geneSymbolmore1 = paste(unique(stringr::str_trim(pos_genes)), collapse = ", "),
      geneSymbolmore2 = ""  # 负向先留空，更稳
    )
  }))
  db_path <- file.path(tempdir(), paste0("ScTypeDB_", tissue, ".xlsx"))
  openxlsx::write.xlsx(db_df, db_path, overwrite = TRUE)
  
  # ---- 构建 Seurat（与 scCATCH 风格对齐；v5 使用 data layer）----
  seu <- Seurat::CreateSeuratObject(counts = expr_mat, assay = "RNA", project = "Simulated")
  SeuratObject::DefaultAssay(seu) <- "RNA"
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = expr_mat)
  
  if (use_all_genes) {
    VariableFeatures(seu) <- rownames(seu)   # 全基因与 scCATCH 的可选项一致
  } else {
    seu <- Seurat::FindVariableFeatures(seu, assay = "RNA", selection.method = "vst",
                                        nfeatures = nfeatures, layer = "data")
  }
  
  seu <- Seurat::ScaleData(seu, features = VariableFeatures(seu))
  seu <- Seurat::RunPCA(seu, features = VariableFeatures(seu))
  seu <- Seurat::FindNeighbors(seu)
  seu <- Seurat::FindClusters(seu)
  # UMAP 是否画图留给你外部控制；这里不跑/不画
  
  # ---- scType 打分（逐细胞）----
  suppressMessages({
    source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
    source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
  })
  gs_list <- gene_sets_prepare(db_path, tissue)
  
  # 兼容 v4/v5：取 scaled 矩阵
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seu[["RNA"]])))
  scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seu[["RNA"]]$scale.data) else as.matrix(seu[["RNA"]]@scale.data)
  
  es.max <- sctype_score(
    scRNAseqData = scRNAseqData_scaled,
    scaled = TRUE,
    gs  = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  # es.max: 行=cell types, 列=细胞
  
  # ---- 逐细胞唯一标签（which.max；可选阈值）----
  type_names <- rownames(es.max)
  max_idx <- apply(es.max, 2, which.max)
  max_val <- vapply(seq_along(max_idx), function(i) es.max[max_idx[i], i], numeric(1))
  labels_cell <- type_names[max_idx]
  names(labels_cell) <- colnames(es.max)
  
  if (!is.null(min_score)) {
    labels_cell[max_val < min_score] <- "Unknown"
  }
  seu[[out_name_cell]] <- labels_cell[colnames(seu)]
  
  # ---- 按 cluster 汇总（每 cluster 一个唯一标签），再扩散到细胞 ----
  agg_fun <- switch(cluster_method,
                    sum    = function(m) rowSums(m),
                    mean   = function(m) rowMeans(m),
                    median = function(m) apply(m, 1, median))
  
  clu <- Idents(seu)
  clusters <- sort(unique(clu))
  per_cluster <- lapply(clusters, function(cl){
    sel <- names(clu)[clu == cl]
    sub <- es.max[, sel, drop = FALSE]
    sc  <- agg_fun(sub)
    idx <- which.max(sc)  # 平局取第一个
    data.frame(
      cluster = as.character(cl),
      type    = names(sc)[idx],
      score   = unname(sc[idx]),
      ncells  = length(sel),
      stringsAsFactors = FALSE
    )
  })
  cluster_table <- do.call(rbind, per_cluster)
  
  lab_map <- setNames(cluster_table$type, cluster_table$cluster)
  labels_cluster <- unname(lab_map[ as.character(clu) ])
  names(labels_cluster) <- names(clu)
  seu[[out_name_clu]] <- labels_cluster[colnames(seu)]
  
  # ---- 组装 truth + 指标，与 scCATCH 风格一致 ----
  if (!is.null(result$label)) {
    truth_all <- paste0("Type", result$label)
    truth <- truth_all[idx]
  } else {
    truth <- rep(NA_character_, length(labels_cell))
    warning("result$label 不存在，无法计算准确度；仅返回 scType 结果。")
  }
  
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
  
  acc_result <-
    if (!is.null(accuracy_fun) && all(!is.na(truth))) {
      accuracy_fun(labels_cell, truth)
    } else if (all(!is.na(truth))) {
      .fallback_metrics(labels_cell, truth)
    } else {
      NA
    }
  
  list(
    seurat              = seu,
    predictions         = labels_cell,     # 逐细胞唯一标签（与你 scCATCH 的 predictions 同名）
    predictions_cluster = labels_cluster,  # 每个cluster唯一标签，扩散回细胞
    truth               = truth,
    accuracy_result     = acc_result,
    score_matrix        = es.max,
    cluster_table       = cluster_table,
    db_df               = db_df,
    db_path             = db_path
  )
}

all_metrics_ScType <- list()
for (i in 1:20) {
  result <- result_list[[i]]
  res <- run_scType_with_accuracy(result, tissue = "Custom", use_all_genes = TRUE)
  
  # 与你 scCATCH 的评估写法一致：
  true_labels <- res$truth
  annotation_labels <- res$predictions  # 逐细胞唯一标签
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  
  all_classes <- union(unique(annotation_labels), unique(true_labels))
  annotation_labels <- factor(annotation_labels, levels = all_classes)
  true_labels <- factor(true_labels, levels = all_classes)
  
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
  
  all_metrics_ScType[[i]] <- list(
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
  
  cat("scType results for group", i, ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
}