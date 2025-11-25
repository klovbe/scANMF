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

Muraro <- MuraroPancreasData()
# 提取 label 列（细胞类型标签）
Muraro_label <- colData(Muraro)$label
# 获取需要去除的索引（NA 和 "unclear"）
indices_to_remove <- which(is.na(Muraro_label) | Muraro_label == "unclear")
# 从表达矩阵中去除这些样本
Muraro_filtered <- Muraro[, -indices_to_remove]
# 从 Muraro_label 中去除这些标签
Muraro_label_filtered <- Muraro_label[-indices_to_remove]


Muraro_filtered <- CreateSeuratObject(counts(Muraro_filtered),project = "Muraro",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
Muraro_filtered[["percent.mt"]] <- PercentageFeatureSet(Muraro_filtered, pattern = "^MT-")
Muraro_filtered <- subset(Muraro_filtered, subset = percent.mt < 20)
# 清理基因名
gene_names <- sub("--chr.*", "", rownames(Muraro_filtered))
Muraro_filtered <- Muraro_filtered[!duplicated(gene_names), ]
rownames(Muraro_filtered) <- gene_names[!duplicated(gene_names)]

Muraro_filtered <- NormalizeData(Muraro_filtered)
Muraro_filtered <- FindVariableFeatures(Muraro_filtered)
Muraro_filtered<- ScaleData(Muraro_filtered, features = rownames(Muraro_filtered))


highly_variable_genes_Muraro <- VariableFeatures(Muraro_filtered, n=1000)
expression_matrix_Muraro <- GetAssayData(Muraro_filtered, assay = "RNA")
# expression_matrix_Muraro <- t(apply(expression_matrix_Muraro, 1, function(x) x / (sd(x) + 1e-10)))
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
union_vec <- toupper(union(highly_variable_genes_Muraro, unique_vector))
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

expression_matrix_Muraro<-as.matrix(expression_matrix_Muraro)

# 高变基因+标记基因
#基因表达矩阵与高变基因+标记基因相交部分
common_rows <- intersect(toupper(rownames(expression_matrix_Muraro)), toupper(rownames(markers)))
# common_rows 现在包含了两个矩阵中行名相同的行的标识符
common_rows1 <- rownames(expression_matrix_Muraro)[toupper(rownames(expression_matrix_Muraro)) %in% common_rows]
expression_matrix_common <- expression_matrix_Muraro[common_rows1, , drop = FALSE]
rownames(expression_matrix_common)<-toupper(rownames(expression_matrix_common)) 
markers_common<- markers[common_rows, ]
colnames(markers_common)
#zero_cols <- which(colSums(markers_common == 0) == nrow(markers_common))
# 删除全为零的列
#markers_common<- markers_common[, -zero_cols]
mapping <- c(
  "acinar" = "Acinar cells",
  "alpha" = "Alpha cells",
  "beta" = "Beta cells",
  "delta" = "Delta cells",
  "duct" = "Ductal cells",
  "epsilon" = "Epsilon cells",
  "pp" = "Gamma (PP) cells",
  "endothelial" = "Endothelial cells",
  "mesenchymal" = "Mesenchymal cells"
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
Muraro_label_filtered <- ifelse(Muraro_label_filtered %in% names(mapping), mapping[Muraro_label_filtered], Muraro_label_filtered)

#——————————————————————————————————
#标签信息
unique_labels <- unique(Muraro_label_filtered)
num_labels <- length(unique_labels)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells <- length(Muraro_label_filtered)
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells) {
  labels <- Muraro_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数


########################
#########以下是sccatch和sctype
seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
result_scCATCH <- list()
result_sctype <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Muraro_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  seu <- Seurat::RunPCA(Muraro_filtered[,!retain_mask], features = VariableFeatures(Muraro_filtered), verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, verbose = FALSE)
  # seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, graph = NULL, nn.name = NULL, features = NULL) 
  # DimPlot(seu, reduction = "umap", label = TRUE)


  obj <- scCATCH::createscCATCH(
    data    = as.matrix(GetAssayData(seu, assay="RNA", layer="data")),
    cluster = as.character(Idents(seu))
  )
  cellmatch_new <- cellmatch[cellmatch$species == "Human" & cellmatch$tissue %in% c("Pancreatic islet","Pancreas","Pancreatic acinar tissue"), ]
  obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
  # obj <- scCATCH::findmarkergene(
  #   object  = obj,
  #   species = "Human",
  #   if_use_custom_marker = TRUE, 
  #   marker= cellmatch_new
  #   # cell_min_pct = 0.02,        # 放宽 marker 表达比例
  #   # logfc = 0.1, pvalue = 0.1     # 放宽显著性阈值
  # )
  
  # 5) 细胞
  obj <- scCATCH::findcelltype(obj)
  # 4) 保护：确认这次有候选
  stopifnot(!is.null(obj@celltype), nrow(obj@celltype) > 0)
  
  
  # 6) 映射逐细胞标签
  ct <- as.data.frame(obj@celltype)
  col_cluster <- intersect(c("cluster","Cluster","clusters"), colnames(ct))[1]
  col_type    <- intersect(c("celltype","cell_type","cell.type","Celltype","CellType"), colnames(ct))[1]
  stopifnot(!is.na(col_cluster), !is.na(col_type))
  cl_map <- setNames(ct[[col_type]], ct[[col_cluster]])
  seu$scCATCH_celltype <- unname(cl_map[as.character(Idents(seu))])
  pred_celltype <- seu$scCATCH_celltype
  
  cat("per-cluster 注释条目：", nrow(ct), "\n",
      "逐细胞 NA 比例：", mean(is.na(seu$scCATCH_celltype)), "\n")
  
  mapping <- c(
    "Acinar Cell" = "Acinar cells",
    "Alpha Cell" = "Alpha cells",
    "Beta Cell" = "Beta cells",
    "Delta Cell" = "Delta cells",
    "PP Cell" = "Gamma (PP) cells"
  )
  
  # 替换 truth 中的名称（不在 mapping 的保持原样）
  pred_aligned <- ifelse(pred_celltype %in% names(mapping), mapping[pred_celltype], pred_celltype)
  pred_aligned[is.na(pred_aligned)] <- "unknown"
  
  ## 去掉名字属性（保持向量干净）
  pred_aligned <- unname(pred_aligned)

  ## 检查结果
  table(pred_aligned, Muraro_label_filtered[!retain_mask])
  mean(pred_aligned == Muraro_label_filtered[!retain_mask])
  result_scCATCH[[i]] <- metric_cal(Muraro_label_filtered[!retain_mask], pred_aligned)


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
  
  
  out_name_cell = "ScType_single"
  out_name_clu  = "ScType_cluster"
  cluster_method = "sum"  # 聚合方式 c("sum","mean","median")

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
  pred_aligned <- labels_cluster
  pred_aligned[is.na(pred_aligned)] <- "unknown"
  
  ## 去掉名字属性（保持向量干净）
  pred_aligned <- unname(pred_aligned)
  
  ## 检查结果
  result_sctype[[i]] <- metric_cal(Muraro_label_filtered[!retain_mask], pred_aligned)
}
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
result_scanmf <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Muraro_label_filtered)
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
  
  # +   alpha<-10000000
  # +   beta<-1000
  # +   lambda<-10000000
  alpha<-10000
  beta<-10
  lambda<-10000
  
  resultMuraro<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
  cell_type_matrix_Muraro  <-resultMuraro[[2]]
  marker_Muraro<-resultMuraro[[1]]
  
  max_col_indices <- apply(cell_type_matrix_Muraro, 1, which.max)
  ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
  normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix_Muraro ), ncol = ncol(cell_type_matrix_Muraro))
  for (j in seq_along(max_col_indices)) {
    normalized_matrix[j, max_col_indices[j]] <- 1
  }
  cell_type_matrix <- normalized_matrix
  colnames(cell_type_matrix)<-colnames(labelMatrix)
  annotation_Muraro<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
  # 计算注释准确率
  true_labels<-Muraro_label_filtered[!retain_mask]
  annotation_labels<-annotation_Muraro[!retain_mask]
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  # 输出准确率
  print(accuracy)
  table(annotation_labels, true_labels)
  result_scanmf[[i]] <- metric_cal(true_labels, annotation_labels)
}


#######################
###singleR
# 设置保留标签信息的比例
result_sigleR <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Muraro_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  x <- expression_matrix_Muraro[,retain_mask]
  y <- expression_matrix_Muraro[,!retain_mask]
  true_labels<-Muraro_label_filtered[!retain_mask]
  pred <- SingleR(y,x,labels = Muraro_label_filtered[retain_mask])
  accuracy <- sum(true_labels == pred$labels) / length(true_labels)
  # 输出准确率
  print(accuracy)
  table(pred$labels, true_labels)  
  result_sigleR[[i]] <-   metric_cal(true_labels, pred$labels)
}



#######################
#####scpred
result_scpred <- list()
for (i in 1:10) {
  # 尝试执行每次采样与注释
  set.seed(seed_vec[i])
  n <- length(Muraro_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  result <- tryCatch({
    x <- GetAssayData(Muraro_filtered, assay = "RNA",layer = "counts")
    ref <- x[,retain_mask]
    query<-x[,!retain_mask]
    ref_label <- Muraro_label_filtered[retain_mask]
    query_label <- Muraro_label_filtered[!retain_mask]
    
    ref <- CreateSeuratObject(ref, project = "Muraro", min.cells = 3, min.features = 200)
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
    ref <- ScaleData(ref, features = rownames(ref), verbose = FALSE)
    ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 10, reduction.key = "pca_", verbose = FALSE, approx = FALSE)
    ref@meta.data$cell_type <- ref_label
    DefaultAssay(ref) <- "RNA"
    
    
    query <- CreateSeuratObject(query, project = "Muraro")
    query <- NormalizeData(query)
    DefaultAssay(query) <- "RNA"
    ref <- getFeatureSpace(ref, "cell_type")
    ref <- trainModel(ref)
    query <- scPredict(query, ref)
    # DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
    pred_scpred <- query@meta.data$scpred_prediction  
    metric_cal(query_label, pred_scpred)
  },
  error = function(e) {
    message("⚠️ 第 ", i, " 次运行出错：", conditionMessage(e))
    return(NA)  # 出错时返回 NA，不中断循环
  })
  
  result_scpred[[i]] <- result
}

Muraro<- list(
  result_scanmf = result_scanmf,
  result_scCATCH = result_scCATCH,
  result_scpred = result_scpred,
  result_sctype = result_sctype,
  result_sigleR = result_sigleR
)

save(Muraro, file = "Muraro_in.RData")
save(result_scanmf, result_scCATCH, result_scpred, result_sctype, result_sigleR, file = "Muraro_in.RData")

