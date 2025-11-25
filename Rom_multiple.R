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
# Load dataset mouse
Romanov <- RomanovBrainData()
unique(colData(Romanov)$'level1 class')
Romanov_label_filtered <- colData(Romanov)$'level1 class'

library(homologene)

# 获取人鼠同源表
orthologs <- homologene::homologeneData
colnames(orthologs)
# 通常列名为: "HomoloGeneID" "Taxonomy_ID" "Gene_ID" "Symbol" "Protein_gi" "Protein_accession"

# 筛选出人类(9606)和小鼠(10090)
orthologs_hm <- subset(orthologs, Taxonomy %in% c(9606, 10090))

# 转为宽表形式（人类 ↔ 小鼠）
orthologs_hm <- reshape(
  orthologs_hm[, c("HID", "Gene.Symbol", "Taxonomy")],
  timevar = "Taxonomy",
  idvar = "HID",
  direction = "wide"
)

# 重命名列
colnames(orthologs_hm) <- c("HomoloGeneID", "human_symbol", "mouse_symbol")

# 查看部分映射
head(orthologs_hm)

# 提取映射关系
mapping_df <- orthologs_hm[, c("mouse_symbol", "human_symbol")]

# 取出小鼠表达矩阵行名
mouse_genes <- rownames(Romanov)

# 映射
human_mapped_genes <- mapping_df$human_symbol[
  match(mouse_genes, mapping_df$mouse_symbol)
]

# 替换行名
rownames(Romanov) <- human_mapped_genes

# 去掉没有匹配的人类基因
Romanov <- Romanov[!is.na(rownames(Romanov)), ]
Romanov <- Romanov[rownames(Romanov) != "", ]
Romanov <- Romanov[!duplicated(rownames(Romanov)), ]

Romanov_filtered <- CreateSeuratObject(counts(Romanov),project = "Romanov",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
Romanov_filtered[["percent.mt"]] <- PercentageFeatureSet(Romanov_filtered, pattern = "^MT-")
Romanov_filtered <- subset(Romanov_filtered, subset = percent.mt < 20)
Romanov_filtered <- NormalizeData(Romanov_filtered)
Romanov_filtered <- FindVariableFeatures(Romanov_filtered)
Romanov_filtered<- ScaleData(Romanov_filtered, features = rownames(Romanov_filtered))
highly_variable_genes_Romanov <- VariableFeatures(Romanov_filtered, n=1000)
expression_matrix_Romanov <- GetAssayData(Romanov_filtered, assay = "RNA")
expression_matrix_Romanov <- as.matrix(expression_matrix_Romanov)

checkGeneSymbols <- function(x, ...) {
  data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
}

source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
# DB file
library(openxlsx)
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare(db_, tissue)
gs = gs_list$gs_positive
gs$`Cancer cells` <- NULL
gs$`Cancer stem cells` <- NULL
gs_list1 = gene_sets_prepare(db_, "Hippocampus")
gs1 = gs_list1$gs_positive
gs_new <- c(gs, gs1[c('Ependymal',"Smooth Muscle")])

# assign cell types
matrix_from_list <- do.call(cbind, gs_new)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
union_vec <- toupper(union(highly_variable_genes_Romanov, unique_vector))
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

# 高变基因+标记基因
#基因表达矩阵与高变基因+标记基因在高边基因内部的部分
common_rows <- intersect(toupper(rownames(expression_matrix_Romanov)), toupper(rownames(markers)))
# common_rows 现在包含了两个矩阵中行名相同的行的标识符
common_rows1 <- rownames(expression_matrix_Romanov)[toupper(rownames(expression_matrix_Romanov)) %in% common_rows]
expression_matrix_common <- expression_matrix_Romanov[common_rows1, , drop = FALSE]
rownames(expression_matrix_common)<-toupper(rownames(expression_matrix_common)) 
markers_common<- markers[common_rows, ]
colnames(markers_common)
#zero_cols <- which(colSums(markers_common == 0) == nrow(markers_common))
# 删除全为零的列
#markers_common<- markers_common[, -zero_cols]
unique(Romanov_label_filtered)

col1 <- markers_common[, 2]
col2 <- markers_common[, 3]
col3 <- markers_common[, 5]
col4 <- markers_common[, 6]
col5 <- markers_common[, 7]
col6 <- markers_common[, 9]
col7 <- markers_common[, 21]

# 创建新的合并列
merged_col <- ifelse(col1 == col2, col2, 1)
merged_col<-ifelse(merged_col == col3, col3, 1)
merged_col <- ifelse(merged_col == col4, col4, 1)
merged_col<-ifelse(merged_col == col5, col5, 1)
merged_col <- ifelse(merged_col == col6, col6, 1)
merged_col<-ifelse(merged_col == col7, col7, 1)
Neurons<-merged_col
# 合并后的新矩阵
# markers_common <-markers_common[, -c(2,3,5,6,7,9,21,23,24)]
markers_common <-markers_common[, -c(2,3,5,6,7,9,21)]
markers_common <- cbind(markers_common, Neurons)


Romanov_label_filtered <- gsub("ependymal", "Ependymal", Romanov_label_filtered)
Romanov_label_filtered <- gsub("vsm", "Smooth Muscle", Romanov_label_filtered)
Romanov_label_filtered <- gsub("neurons", "Neurons", Romanov_label_filtered)
Romanov_label_filtered <- gsub("oligos", "Oligodendrocytes", Romanov_label_filtered)
Romanov_label_filtered <- gsub("microglia", "Microglial cells", Romanov_label_filtered)
Romanov_label_filtered <- gsub("endothelial", "Endothelial cells", Romanov_label_filtered)
Romanov_label_filtered <- gsub("astrocytes", "Astrocytes", Romanov_label_filtered)


#——————————————————————————————————
#标签信息
unique_labels <- unique(Romanov_label_filtered)
num_labels <- length(unique_labels)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells <- length(Romanov_label_filtered)
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))
# 填充标签矩阵
for (i in 1:num_cells) {
  labels <- Romanov_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数


########################
#########以下是sccatch和sctype
# 设置保留比例
# 设置保留标签信息的比例
seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
result_scCATCH <- list()
result_sctype <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Romanov_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  # 应用掩码，保留的细胞标签不变，其余细胞整行变为 0
  seu <- Seurat::RunPCA(Romanov_filtered[,!retain_mask], features = VariableFeatures(Romanov_filtered), verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, verbose = FALSE)
  seu_data <- GetAssayData(seu, assay="RNA", layer="data")
  # seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, graph = NULL, nn.name = NULL, features = NULL) 
  # DimPlot(seu, reduction = "umap", label = TRUE)
  
  # 取出小鼠表达矩阵行名
  human_genes <- rownames(seu_data)
  
  # 映射
  mouse_mapped_genes <- mapping_df$mouse_symbol[
    match(human_genes, mapping_df$human_symbol)
  ]
  
  # 替换行名
  rownames(seu_data) <- mouse_mapped_genes
  
  
  
  obj <- scCATCH::createscCATCH(
    data    = as.matrix(seu_data),
    cluster = as.character(Idents(seu))
  )

  ## 严格用系统库分支（不要传 marker；只传 species/tissue）
  # 尝试 “Pancreatic islet”（更贴内分泌）
  # cellmatch_new <- cellmatch[cellmatch$species == "Mouse" & cellmatch$tissue %in% c("Brain"), ]
  # obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
  obj <- findmarkergene(object = obj, species = "Mouse", marker = cellmatch, tissue = "Brain")
  
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
  map_unique <- data.frame(
    pred = c(
      "Oligodendrocyte",
      "Oligodendrocyte, spiral ganglion neuron",
      "Astrocyte",
      "Microglial Cell",
      "Endothelial Cell",
      "Type IC spiral ganglion neuron",
      "spiral ganglion neuron",
      "Type I spiral ganglion neuron",
      "Type II spiral ganglion neuron",
      "Ependymal Cell" 
    ),
    mapped_truth = c(
      "Oligodendrocytes",
      "Oligodendrocytes",
      "Astrocytes",
      "Microglia cells",
      "Endothelial cells",
      "Neurons",
      "Neurons",
      "Neurons",
      "Neurons",
      "Ependymal"
    ),
    stringsAsFactors = FALSE
  )
  
  
  ## 生成映射字典
  map_vec <- setNames(map_unique$mapped_truth, map_unique$pred)
  
  
  ## 用 ifelse 判断，能对得上就改，否则保留原值
  pred_aligned <- ifelse(pred_celltype %in% names(map_vec),
                         map_vec[pred_celltype],
                         pred_celltype)
  pred_aligned[is.na(pred_aligned)] <- "unknown"
  
  ## 去掉名字属性（保持向量干净）
  pred_aligned <- unname(pred_aligned)

  table(pred_aligned, Romanov_label_filtered[!retain_mask])
  mean(pred_aligned == Romanov_label_filtered[!retain_mask])
  result_scCATCH[[i]] <- metric_cal(Romanov_label_filtered[!retain_mask], pred_aligned)


  #sctype
  checkGeneSymbols <- function(x, ...) {
    data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
  }
  # ---- scType 打分（逐细胞）----
  suppressMessages({
    source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
    source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
  }) 
  # DB file
  db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
  tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
  
  # prepare gene sets
  gs_list <- gene_sets_prepare(db_, tissue)
  
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
  
  map_unique <- data.frame(
    pred = c(
      "GABAergic neurons",
      "Glutamatergic neurons",
      "Immature neurons",
      "Mature neurons"
    ),
    mapped_truth = c(
      "Neurons",
      "Neurons",
      "Neurons",
      "Neurons"
    ),
    stringsAsFactors = FALSE
  )
  
  
  ## 生成映射字典
  map_vec <- setNames(map_unique$mapped_truth, map_unique$pred)
  
  
  ## 用 ifelse 判断，能对得上就改，否则保留原值
  pred_aligned <- ifelse(labels_cluster %in% names(map_vec),
                         map_vec[labels_cluster],
                         labels_cluster)
  pred_aligned[is.na(pred_aligned)] <- "unknown"
  
  ## 去掉名字属性（保持向量干净）
  pred_aligned <- unname(pred_aligned)
  
  ## 检查结果
  result_sctype[[i]] <- metric_cal(Romanov_label_filtered[!retain_mask], pred_aligned)
}

###############################################
#mutiple rounds
########################################################
#Laplace
# 转置：列为样本，行为特征
expression_matrix_common <- as.matrix(expression_matrix_common)
X_t <- t(expression_matrix_common)
#KNN+Gauss
k <- 1000
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




# 设置保留比例
# 设置保留标签信息的比例
seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
result_scanmf <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Romanov_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  # 应用掩码，保留的细胞标签不变，其余细胞整行变为 0
  masked_label_matrix <- label_matrix * retain_mask  # 逐行运算，FALSE 位置的整行变为 0
  colnames(masked_label_matrix)<-unique_labels
  colnames(label_matrix)<-colnames(masked_label_matrix)
  #标签和标记基因对应的细胞类型的并集
  unique_labellist <- unique(c(colnames(markers_common),colnames(label_matrix)))
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
  
  # +   alpha<-500
  # +   beta<-10
  # +   lambda<-1000000
  alpha<-10000
  beta<-10
  lambda<-10000
  
  resultRomanov<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
  cell_type_matrix_Romanov  <-resultRomanov[[2]]
  marker_Romanov<-resultRomanov[[1]]
  
  max_col_indices <- apply(cell_type_matrix_Romanov, 1, which.max)
  ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
  normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix_Romanov ), ncol = ncol(cell_type_matrix_Romanov))
  for (j in seq_along(max_col_indices)) {
    normalized_matrix[j, max_col_indices[j]] <- 1
  }
  cell_type_matrix <- normalized_matrix
  colnames(cell_type_matrix)<-colnames(labelMatrix)
  annotation_Romanov<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
  # 计算注释准确率
  true_labels<-Romanov_label_filtered[!retain_mask]
  annotation_labels<-annotation_Romanov[!retain_mask]
  accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
  # 输出准确率
  print(accuracy)
  table(annotation_labels, true_labels)
  result_scanmf[[i]] <- metric_cal(true_labels, annotation_labels)
}

#######################
###singleR
# 设置保留标签信息的比例
retain_ratio <- 0.02
result_sigleR <- list()
for (i in 1:10){
  set.seed(seed_vec[i])
  n <- length(Romanov_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  x <- expression_matrix_Romanov[,retain_mask]
  y <- expression_matrix_Romanov[,!retain_mask]
  true_labels<-Romanov_label_filtered[!retain_mask]
  pred <- SingleR(y,x,labels = Romanov_label_filtered[retain_mask])
  accuracy <- sum(true_labels == pred$labels) / length(true_labels)
  # 输出准确率
  print(accuracy)
  table(pred$labels, true_labels)  
  result_sigleR[[i]] <-   metric_cal(true_labels, pred$labels)
}



#######################
#####scpred

#######################
#####scpred
result_scpred <- list()
for (i in 1:10) {
  # 尝试执行每次采样与注释
  set.seed(seed_vec[i])
  n <- length(Romanov_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  result <- tryCatch({
    x <- GetAssayData(Romanov_filtered, assay = "RNA",layer = "counts")
    ref <- x[,retain_mask]
    query<-x[,!retain_mask]
    ref_label <- Romanov_label_filtered[retain_mask]
    query_label <- Romanov_label_filtered[!retain_mask]
    
    ref <- CreateSeuratObject(ref, project = "Romanov", min.cells = 3, min.features = 200)
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
    ref <- ScaleData(ref, features = rownames(ref), verbose = FALSE)
    ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 10, reduction.key = "pca_", verbose = FALSE, approx = FALSE)
    ref@meta.data$cell_type <- ref_label
    DefaultAssay(ref) <- "RNA"
    
    
    query <- CreateSeuratObject(query, project = "Romanov")
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

Romanov <- list(
  result_scanmf = result_scanmf,
  result_scCATCH = result_scCATCH,
  result_scpred = result_scpred,
  result_sctype = result_sctype,
  result_sigleR = result_sigleR
)



save(Romanov,file = "Romanov_in.RData")

load("Romanov_in.RData")


