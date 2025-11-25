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

################################################################
#cross dataset
# Load dataset mouse
human <- BaronPancreasData('human')
mouse <- BaronPancreasData('mouse')
human_label <- colData(human)$label
mouse_label <- colData(mouse)$label
mouse <- genem2h(mouse)

# 去掉没有匹配的人类基因
mouse <- mouse[!is.na(rownames(mouse)), ]
mouse <- mouse[rownames(mouse) != "", ]
mouse <- mouse[!duplicated(rownames(mouse)), ]


# 找到 fetal_quiescent 和 hybrid 的索引 
indices_to_remove <- which(mouse_label %in% c("B_cell", "immune_other"))
# indices_to_remove <- which(mouse_label %in% c("B_cell", "immune_other", "macrophage","t_cell","T_cell"))
# indices_to_remove <- which(mouse_label == "hybrid")
# 删除对应的标签信息
mouse_label_filtered <- mouse_label[-indices_to_remove]
mouse_filtered <- mouse[, -indices_to_remove]
indices_to_remove <- which(human_label %in% c("B_cell", "immune_other"))
# indices_to_remove <- which(human_label %in% c("B_cell", "immune_other", "macrophage","t_cell","T_cell"))
human_label_filtered <- human_label
human_filtered <- human

mouse_filtered <- CreateSeuratObject(counts(mouse_filtered),project = "mouse",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
mouse_filtered[["percent.mt"]] <- PercentageFeatureSet(mouse_filtered, pattern = "^MT-")
mouse_filtered <- subset(mouse_filtered, subset = percent.mt < 20)
mouse_filtered <- NormalizeData(mouse_filtered)
mouse_filtered <- FindVariableFeatures(mouse_filtered)
mouse_filtered<- ScaleData(mouse_filtered, features = rownames(mouse_filtered))
highly_variable_genes_mouse <- VariableFeatures(mouse_filtered, n=1000)
expression_matrix_mouse <- GetAssayData(mouse_filtered, assay = "RNA")
expression_matrix_mouse<-as.matrix(expression_matrix_mouse)



human_filtered <- CreateSeuratObject(counts(human_filtered),project = "human",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
human_filtered[["percent.mt"]] <- PercentageFeatureSet(human_filtered, pattern = "^MT-")
human_filtered <- subset(human_filtered, subset = percent.mt < 20)
human_filtered <- NormalizeData(human_filtered)
human_filtered <- FindVariableFeatures(human_filtered)
human_filtered<- ScaleData(human_filtered, features = rownames(human_filtered))
highly_variable_genes_human <- VariableFeatures(human_filtered, n=1000)
expression_matrix_human <- GetAssayData(human_filtered, assay = "RNA")

mapping <- c(
  "acinar" = "Acinar cells",
  "alpha" = "Alpha cells",
  "beta" = "Beta cells",
  "delta" = "Delta cells",
  "ductal" = "Ductal cells",
  "endothelial"  = "Endothelial cells",
  "epsilon" = "Epsilon cells",
  "gamma" = "Gamma (PP) cells",
  "PSC" = "Pancreatic stellate cells",
  "mast" = "Mast cells",
  "schwann" = "Peri-islet Schwann cells",
  "activated_stellate" = "Pancreatic stellate cells",
  "quiescent_stellate" = "Pancreatic stellate cells",
  "macrophage" = "Immune system cells",
  "t_cell" ="Immune system cells",
  "T_cell"="Immune system cells"
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
human_label_filtered <- ifelse(human_label_filtered %in% names(mapping), mapping[human_label_filtered], human_label_filtered)
mouse_label_filtered <- ifelse(mouse_label_filtered %in% names(mapping), mapping[mouse_label_filtered], mouse_label_filtered)



checkGeneSymbols <- function(x, ...) {
  data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
}

source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
# DB file
library(openxlsx)
tissue = "Pancreas" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare(db_, tissue)
gs = gs_list$gs_positive
gs$`Cancer cells` <- NULL
gs$`Cancer stem cells` <- NULL

z <- unique(intersect(union(toupper(highly_variable_genes_human),toupper(highly_variable_genes_mouse)),intersect(toupper(rownames(expression_matrix_human)), toupper(rownames(expression_matrix_mouse)))))
# assign cell types
matrix_from_list <- do.call(cbind, gs)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
unique_vector <- unique(intersect(unique_vector,intersect(toupper(rownames(expression_matrix_human)), toupper(rownames(expression_matrix_mouse)))))
union_vec <- toupper(union(z, unique_vector))
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
rownames(expression_matrix_mouse) <- toupper(rownames(expression_matrix_mouse))
rownames(expression_matrix_human) <- toupper(rownames(expression_matrix_human))
expression_matrix_common_mouse <- expression_matrix_mouse[union_vec, , drop = FALSE]
expression_matrix_common_human <- expression_matrix_human[union_vec, , drop = FALSE]



expression_matrix_common <- cbind(expression_matrix_common_mouse, expression_matrix_common_human)
markers_common<- markers
num_cells_mouse <- length(mouse_label_filtered)
num_cells_human <- length(human_label_filtered)
num_cells <- num_cells_human+num_cells_mouse
#——————————————————————————————————
#标签信息
unique_labels <- unique(mouse_label_filtered)
num_labels <- length(unique_labels)
# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells_mouse) {
  labels <- mouse_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数# 设置保留比例

scanmf_timesh <- c()

for (kk in 1:1) {
  t0 <- Sys.time()

# 查看结果
colnames(label_matrix)<-unique_labels
#标签和标记基因对应的细胞类型的并集
unique_labellist <- unique(c(colnames(markers_common),colnames(label_matrix)))
#标签以外的细胞类型
a<- setdiff(unique_labellist, colnames(label_matrix))
#额外的细胞类型矩阵
label_matrix1 <- matrix(0, nrow = num_cells, ncol =length(a) )
colnames(label_matrix1)<-a
#总的含有标签信息细胞类型矩阵
labelMatrix<-cbind(label_matrix1 ,label_matrix)
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

#Laplace
expression_matrix_common<-as.matrix(expression_matrix_common)
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
# 若需对称归一化，可用：
# D_inv_sqrt <- diag(1 / sqrt(diag(D) + 1e-10))
# L <- diag(n) - D_inv_sqrt %*% Q %*% D_inv_sqrt


#0.933333 alpha/40
alpha<-1000
beta<-10
lambda<-1000

result<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
cell_type_matrix  <-result[[2]]
marker<-result[[1]]

max_col_indices <- apply(cell_type_matrix, 1, which.max)
##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix), ncol = ncol(cell_type_matrix))
for (j in seq_along(max_col_indices)) {
  normalized_matrix[j, max_col_indices[j]] <- 1
}
cell_type_matrix <- normalized_matrix
colnames(cell_type_matrix)<-colnames(labelMatrix)
annotation_labels<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])

t1 <- Sys.time()
scanmf_timesh[kk] <- as.numeric(t1 - t0, units = "secs")
}



######################
#SingleR 0.8982456
singler_timesh <- c()

for (i in 1:1) {
  t0 <- Sys.time()
  
library(SingleR)
pred<- SingleR(test=expression_matrix_human, ref=expression_matrix_mouse, labels=mouse_label_filtered)
pred_singler <- pred$labels
t1 <- Sys.time()
singler_timesh[i] <- as.numeric(t1 - t0, units = "secs")
}



#####scpred 0.3929825
# Load dataset mouse
scpred_timesh <- c()

for (i in 1:1) {
  t0 <- Sys.time()
ref <- GetAssayData(mouse_filtered, assay = "RNA",layer = "counts")
query<-GetAssayData(human_filtered, assay = "RNA",layer = "counts")
ref_label <- mouse_label_filtered
query_label <- human_label_filtered

ref <- CreateSeuratObject(ref, project = "ref", min.cells = 3, min.features = 200)
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref, nfeatures = 2000)
ref <- ScaleData(ref, features = rownames(ref))
ref <- RunPCA(ref, verbose = FALSE)
ref@meta.data$cell_type <- ref_label
DefaultAssay(ref) <- "RNA"


query <- CreateSeuratObject(query, project = "query")
query <- NormalizeData(query)
DefaultAssay(query) <- "RNA"
ref <- getFeatureSpace(ref, "cell_type")
ref <- trainModel(ref)
query <- scPredict(query, ref)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
pred_scpred <- query@meta.data$scpred_prediction  
t1 <- Sys.time()
scpred_timesh[i] <- as.numeric(t1 - t0, units = "secs")
}



#####################
#以human预测mouse
############################
singler_timesm <- c()

for (i in 1:1) {
  t0 <- Sys.time()
  library(SingleR)
  pred<- SingleR(test=expression_matrix_mouse, ref=expression_matrix_human, labels=human_label_filtered)
  pred_singler <- pred$labels
  t1 <- Sys.time()
  singler_timesm[i] <- as.numeric(t1 - t0, units = "secs")
}

#####scpred 0.3929825
# Load dataset mouse
scpred_timesm <- c()

for (i in 1:1) {
  t0 <- Sys.time()
  ref <- GetAssayData(human_filtered, assay = "RNA",layer = "counts")
  query<-GetAssayData(mouse_filtered, assay = "RNA",layer = "counts")
  ref_label <- human_label_filtered
  query_label <- mouse_label_filtered
  
  ref <- CreateSeuratObject(ref, project = "ref", min.cells = 3, min.features = 200)
  ref <- NormalizeData(ref)
  ref <- FindVariableFeatures(ref, nfeatures = 2000)
  ref <- ScaleData(ref, features = rownames(ref))
  ref <- RunPCA(ref, verbose = FALSE)
  ref@meta.data$cell_type <- ref_label
  DefaultAssay(ref) <- "RNA"
  
  
  query <- CreateSeuratObject(query, project = "query")
  query <- NormalizeData(query)
  DefaultAssay(query) <- "RNA"
  ref <- getFeatureSpace(ref, "cell_type")
  ref <- trainModel(ref)
  query <- scPredict(query, ref)
  DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
  pred_scpred <- query@meta.data$scpred_prediction  
  t1 <- Sys.time()
  scpred_timesm[i] <- as.numeric(t1 - t0, units = "secs")
}



#-----------------------------
#scANMF
#——————————————————————————————————
#标签信息
unique_labels <- unique(human_label_filtered)
num_labels <- length(unique_labels)
expression_matrix_common <- cbind(expression_matrix_common_human, expression_matrix_common_mouse)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells_human) {
  labels <- human_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数# 设置保留比例

scanmf_timesm <- c()

for (kk in 1:1) {
  t0 <- Sys.time()
  
  # 查看结果
  colnames(label_matrix)<-unique_labels
  #标签和标记基因对应的细胞类型的并集
  unique_labellist <- unique(c(colnames(markers_common),colnames(label_matrix)))
  #标签以外的细胞类型
  a<- setdiff(unique_labellist, colnames(label_matrix))
  #额外的细胞类型矩阵
  label_matrix1 <- matrix(0, nrow = num_cells, ncol =length(a) )
  colnames(label_matrix1)<-a
  #总的含有标签信息细胞类型矩阵
  labelMatrix<-cbind(label_matrix1 ,label_matrix)
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
  
  #Laplace
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
  # 若需对称归一化，可用：
  # D_inv_sqrt <- diag(1 / sqrt(diag(D) + 1e-10))
  # L <- diag(n) - D_inv_sqrt %*% Q %*% D_inv_sqrt
  
  
  alpha<-100000
  beta<-10
  lambda<-100
  
  result<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
  cell_type_matrix  <-result[[2]]
  marker<-result[[1]]
  
  max_col_indices <- apply(cell_type_matrix, 1, which.max)
  ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
  normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix), ncol = ncol(cell_type_matrix))
  for (j in seq_along(max_col_indices)) {
    normalized_matrix[j, max_col_indices[j]] <- 1
  }
  cell_type_matrix <- normalized_matrix
  colnames(cell_type_matrix)<-colnames(labelMatrix)
  annotation_labels<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
  t1 <- Sys.time()
  scanmf_timesm[kk] <- as.numeric(t1 - t0, units = "secs")
}




#---------------
##sccatch
sccatch_timesh <- c()

for (i in 1:1) {
  t0 <- Sys.time()
human_filtered <- Seurat::RunPCA(human_filtered, features = VariableFeatures(human_filtered), verbose = FALSE)
human_filtered <- Seurat::FindNeighbors(human_filtered, verbose = FALSE)
human_filtered <- Seurat::FindClusters(human_filtered, verbose = FALSE)
obj <- scCATCH::createscCATCH(
  data    = as.matrix(GetAssayData(human_filtered, assay="RNA", layer="data")),
  cluster = as.character(Idents(human_filtered))
)

## 严格用系统库分支（不要传 marker；只传 species/tissue）
# 尝试 “Pancreatic islet”（更贴内分泌）
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
human_filtered$scCATCH_celltype <- unname(cl_map[as.character(Idents(human_filtered))])
pred_celltype <- human_filtered$scCATCH_celltype

cat("per-cluster 注释条目：", nrow(ct), "\n",
    "逐细胞 NA 比例：", mean(is.na(human_filtered$scCATCH_celltype)), "\n")


truth_aligned <- human_label_filtered
mapping <- c(
  "Acinar Cell" = "Acinar cells",
  "Alpha Cell" = "Alpha cells",
  "Beta Cell" = "Beta cells",
  "Delta Cell" = "Delta cells",
  "Epithelial Cell" ="Epithelial cells" 
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
pred_aligned <- ifelse(pred_celltype %in% names(mapping), mapping[pred_celltype], pred_celltype)
pred_aligned[is.na(pred_aligned)] <- "unknown"

## 去掉名字属性（保持向量干净）
pred_aligned <- unname(pred_aligned)
t1 <- Sys.time()
sccatch_timesh[i] <- as.numeric(t1 - t0, units = "secs")
}

sccatch_timesm <- c()

for (i in 1:1) {
  t0 <- Sys.time()

mouse_filtered <- Seurat::RunPCA(mouse_filtered, features = VariableFeatures(mouse_filtered), verbose = FALSE)
mouse_filtered <- Seurat::FindNeighbors(mouse_filtered, verbose = FALSE)
mouse_filtered <- Seurat::FindClusters(mouse_filtered, verbose = FALSE)
obj <- scCATCH::createscCATCH(
  data    = as.matrix(GetAssayData(mouse_filtered, assay="RNA", layer="data")),
  cluster = as.character(Idents(mouse_filtered))
)

## 严格用系统库分支（不要传 marker；只传 species/tissue）
# 尝试 “Pancreatic islet”（更贴内分泌）
cellmatch_new <- cellmatch[cellmatch$species == "Human" & cellmatch$tissue %in% c("Pancreatic islet","Pancreas","Pancreatic acinar tissue"), ]
obj <- findmarkergene(object = obj, if_use_custom_marker = TRUE, marker = cellmatch_new)
# obj <- scCATCH::findmarkergene(
#   object  = obj,
#   species = "mouse",
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
mouse_filtered$scCATCH_celltype <- unname(cl_map[as.character(Idents(mouse_filtered))])
pred_celltype <- mouse_filtered$scCATCH_celltype

cat("per-cluster 注释条目：", nrow(ct), "\n",
    "逐细胞 NA 比例：", mean(is.na(mouse_filtered$scCATCH_celltype)), "\n")


truth_aligned <- mouse_label_filtered
mapping <- c(
  "Acinar Cell" = "Acinar cells",
  "Alpha Cell" = "Alpha cells",
  "Beta Cell" = "Beta cells",
  "Delta Cell" = "Delta cells",
  "Epithelial Cell" ="Epithelial cells" 
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
pred_aligned <- ifelse(pred_celltype %in% names(mapping), mapping[pred_celltype], pred_celltype)
pred_aligned[is.na(pred_aligned)] <- "unknown"

## 去掉名字属性（保持向量干净）
pred_aligned <- unname(pred_aligned)
t1 <- Sys.time()
sccatch_timesm[i] <- as.numeric(t1 - t0, units = "secs")
}

#sctype
sctype_timesm <- c()

for (i in 1:1) {
  t0 <- Sys.time()
out_name_cell = "ScType_single"
out_name_clu  = "ScType_cluster"
min_score     = NULL    # 可选：逐细胞最大分的阈值，小于阈值标 "Unknown"
use_all_genes = TRUE     # TRUE：全部基因设为 VariableFeatures；与 scCATCH 参数一致
nfeatures     = 2000     # use_all_genes=FALSE 时使用
cluster_method = "sum" # 聚合方式 c("sum","mean","median")
# ---- scType 打分（逐细胞）----
suppressMessages({
  source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
  source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
})
# DB file
db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
tissue <- "Pancreas" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# 兼容 v4/v5：取 scaled 矩阵
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(mouse_filtered[["RNA"]])))
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(mouse_filtered[["RNA"]]$scale.data) else as.matrix(mouse_filtered[["RNA"]]@scale.data)

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

min_score <- NULL
if (!is.null(min_score)) {
  labels_cell[max_val < min_score] <- "Unknown"
}
mouse_filtered[[out_name_cell]] <- labels_cell[colnames(mouse_filtered)]

# ---- 按 cluster 汇总（每 cluster 一个唯一标签），再扩散到细胞 ----
agg_fun <- switch(cluster_method,
                  sum    = function(m) rowSums(m),
                  mean   = function(m) rowMeans(m),
                  median = function(m) apply(m, 1, median))

clu <- Idents(mouse_filtered)
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
mouse_filtered[[out_name_clu]] <- labels_cluster[colnames(mouse_filtered)]
truth_aligned <- mouse_label_filtered
pred_aligned <- labels_cluster
pred_aligned[is.na(pred_aligned)] <- "unknown"

## 去掉名字属性（保持向量干净）
pred_aligned <- unname(pred_aligned)
t1 <- Sys.time()
sctype_timesm[i] <- as.numeric(t1 - t0, units = "secs")
}

sctype_timesh <- c()

for (i in 1:1) {
  t0 <- Sys.time()

# 兼容 v4/v5：取 scaled 矩阵
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(human_filtered[["RNA"]])))
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(human_filtered[["RNA"]]$scale.data) else as.matrix(human_filtered[["RNA"]]@scale.data)

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

min_score <- NULL
if (!is.null(min_score)) {
  labels_cell[max_val < min_score] <- "Unknown"
}
human_filtered[[out_name_cell]] <- labels_cell[colnames(human_filtered)]

# ---- 按 cluster 汇总（每 cluster 一个唯一标签），再扩散到细胞 ----
agg_fun <- switch(cluster_method,
                  sum    = function(m) rowSums(m),
                  mean   = function(m) rowMeans(m),
                  median = function(m) apply(m, 1, median))

clu <- Idents(human_filtered)
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
human_filtered[[out_name_clu]] <- labels_cluster[colnames(human_filtered)]
truth_aligned <- human_label_filtered
pred_aligned <- labels_cluster
pred_aligned[is.na(pred_aligned)] <- "unknown"

## 去掉名字属性（保持向量干净）
pred_aligned <- unname(pred_aligned)
t1 <- Sys.time()
sctype_timesh[i] <- as.numeric(t1 - t0, units = "secs")
}

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
Zeisel<- ZeiselBrainData()
unique(colData(Zeisel)$level1class)
Zeisel_label <- colData(Zeisel)$level1class

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
mouse_genes <- rownames(Zeisel)

# 映射
human_mapped_genes <- mapping_df$human_symbol[
  match(mouse_genes, mapping_df$mouse_symbol)
]

# 替换行名
rownames(Zeisel) <- human_mapped_genes

# 去掉没有匹配的人类基因
Zeisel <- Zeisel[!is.na(rownames(Zeisel)), ]
Zeisel <- Zeisel[rownames(Zeisel) != "", ]
Zeisel <- Zeisel[!duplicated(rownames(Zeisel)), ]
Zeisel_label_filtered <- Zeisel_label
Zeisel_filtered <- Zeisel

Zeisel_filtered <- CreateSeuratObject(counts(Zeisel_filtered),project = "Zeisel",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
Zeisel_filtered[["percent.mt"]] <- PercentageFeatureSet(Zeisel_filtered, pattern = "^MT-")
Zeisel_filtered <- subset(Zeisel_filtered, subset = percent.mt < 20)
Zeisel_filtered <- NormalizeData(Zeisel_filtered)
Zeisel_filtered <- FindVariableFeatures(Zeisel_filtered)
Zeisel_filtered<- ScaleData(Zeisel_filtered, features = rownames(Zeisel_filtered))
highly_variable_genes_zeisel <- VariableFeatures(Zeisel_filtered, n=1000)
expression_matrix_Zeisel <- GetAssayData(Zeisel_filtered, assay = "RNA")
expression_matrix_Zeisel<-as.matrix(expression_matrix_Zeisel)


load("Brain_marker.RData")
gs$`Cancer cells` <- NULL
gs$`Cancer stem cells` <- NULL
# assign cell types
matrix_from_list <- do.call(cbind, gs)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
union_vec <- toupper(union(highly_variable_genes_zeisel, unique_vector))
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
common_rows <- intersect(toupper(rownames(expression_matrix_Zeisel)), toupper(rownames(markers)))
# common_rows 现在包含了两个矩阵中行名相同的行的标识符
common_rows1 <- rownames(expression_matrix_Zeisel)[toupper(rownames(expression_matrix_Zeisel)) %in% common_rows]
expression_matrix_common <- expression_matrix_Zeisel[common_rows1, , drop = FALSE]
rownames(expression_matrix_common)<-toupper(rownames(expression_matrix_common)) 
markers_common<- markers[common_rows, ]
colnames(markers_common)
#zero_cols <- which(colSums(markers_common == 0) == nrow(markers_common))
# 删除全为零的列
#markers_common<- markers_common[, -zero_cols]
unique(Zeisel_label_filtered)
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
markers_common <-markers_common[, -c(2,3,5,6,7,9,21)]
markers_common <- cbind(markers_common, Neurons)

# Zeisel_label_filtered <- gsub("interneurons", "GABAergic neurons", Zeisel_label_filtered)
# Zeisel_label_filtered <- gsub("pyramidal SS", "Glutamatergic neurons", Zeisel_label_filtered)
# Zeisel_label_filtered <- gsub("pyramidal CA1", "Glutamatergic neurons", Zeisel_label_filtered)

Zeisel_label_filtered <- gsub("interneurons", "Neurons", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("pyramidal SS", "Neurons", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("pyramidal CA1", "Neurons", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("oligodendrocytes", "Oligodendrocytes", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("microglia", "Microglial cells", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("endothelial-mural", "Endothelial cells", Zeisel_label_filtered)
Zeisel_label_filtered <- gsub("astrocytes_ependymal", "Astrocytes", Zeisel_label_filtered)
#——————————————————————————————————
#标签信息
unique_labels <- unique(Zeisel_label_filtered)
num_labels <- length(unique_labels)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells <- length(Zeisel_label_filtered)
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells) {
  labels <- Zeisel_label_filtered[i]
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
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
sccatch_timesz <- c()
sctype_timesz <- c()
for (i in 1:1){
  t0 <- Sys.time()
  set.seed(seed_vec[i])
  n <- length(Zeisel_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  seu <- Seurat::RunPCA(Zeisel_filtered[,!retain_mask], features = VariableFeatures(Zeisel_filtered), verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, verbose = FALSE)
  seu_data <- GetAssayData(seu, assay="RNA", layer="data")
  # seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, graph = NULL, nn.name = NULL, features = NULL) 
  # DimPlot(seu, reduction = "umap", label = TRUE)
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
      "Astrocyte",
      "Microglial Cell",
      "Endothelial Cell",
      "Type IC spiral ganglion neuron",
      "spiral ganglion neuron"
    ),
    mapped_truth = c(
      "Oligodendrocytes",
      "Astrocytes",
      "Microglia cells",
      "Endothelial cells",
      "Neurons",
      "Neurons"
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
  t1 <- Sys.time()
  sccatch_timesz[i] <- as.numeric(t1 - t0, units = "secs")
  
  
  t0 <- Sys.time()
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
  t1 <- Sys.time()
  sctype_timesz[i] <- as.numeric(t1 - t0, units = "secs")
}

###############################################
#mutiple rounds
########################################################
#Laplace
# 转置：列为样本，行为特征
# 设置保留比例
# 设置保留标签信息的比例
seed_vec <- c(1,11,123,1234,12345,54321,4321,3211,21,0)
# 设置保留标签信息的比例
retain_ratio <- 0.02  # 例如只保留 50% 细胞的标签
scanmf_timesz <- c()
i=1
for (kk in 1:1){
  t0 <- Sys.time()
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
  
  
  
  
  set.seed(seed_vec[1])
  n <- length(Zeisel_label_filtered)
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
  
  # alpha<-500
  # beta<-80
  # lambda<-1000000
  
  alpha<-10000
  beta<-10
  lambda<-10000
  
  resultZeisel<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
  cell_type_matrix_Zeisel  <-resultZeisel[[2]]
  marker_Zeisel<-resultZeisel[[1]]
  
  max_col_indices <- apply(cell_type_matrix_Zeisel, 1, which.max)
  ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
  normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix_Zeisel ), ncol = ncol(cell_type_matrix_Zeisel))
  for (j in seq_along(max_col_indices)) {
    normalized_matrix[j, max_col_indices[j]] <- 1
  }
  cell_type_matrix <- normalized_matrix
  colnames(cell_type_matrix)<-colnames(labelMatrix)
  annotation_Zeisel<-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
  t1 <- Sys.time()
  scanmf_timesz[kk] <- as.numeric(t1 - t0, units = "secs")
}

#######################
###singleR
# 设置保留标签信息的比例
retain_ratio <- 0.02
singler_timesz <- c()
for (i in 1:1){
  t0 <- Sys.time()
  set.seed(seed_vec[i])
  n <- length(Zeisel_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  x <- expression_matrix_Zeisel[,retain_mask]
  y <- expression_matrix_Zeisel[,!retain_mask]
  true_labels<-Zeisel_label_filtered[!retain_mask]
  pred <- SingleR(y,x,labels = Zeisel_label_filtered[retain_mask])
  t1 <- Sys.time()
  singler_timesz[i] <- as.numeric(t1 - t0, units = "secs")
}



#######################
#####scpred
scpred_timesz <- c()
for (i in 1:1) {
  t0 <- Sys.time()
  # 尝试执行每次采样与注释
  set.seed(seed_vec[3])
  n <- length(Zeisel_label_filtered)
  retain_num <- round(n * retain_ratio)
  retain_mask <- rep(FALSE, n)
  retain_mask[sample(seq_len(n), retain_num)] <- TRUE
  result <- tryCatch({
    x <- GetAssayData(Zeisel_filtered, assay = "RNA",layer = "counts")
    ref <- x[,retain_mask]
    query<-x[,!retain_mask]
    ref_label <- Zeisel_label_filtered[retain_mask]
    query_label <- Zeisel_label_filtered[!retain_mask]
    
    ref <- CreateSeuratObject(ref, project = "Zeisel", min.cells = 3, min.features = 200)
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
    ref <- ScaleData(ref, features = rownames(ref), verbose = FALSE)
    ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 10, reduction.key = "pca_", verbose = FALSE, approx = FALSE)
    ref@meta.data$cell_type <- ref_label
    DefaultAssay(ref) <- "RNA"
    
    
    query <- CreateSeuratObject(query, project = "Zeisel")
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
  
  t1 <- Sys.time()
  scpred_timesz[i] <- as.numeric(t1 - t0, units = "secs")
}



# ---- 方法与数据集 ----
methods  <- c("scanmf", "scpred", "sctype", "singler", "sccatch")
datasets <- c("m", "h", "z")   # m=Baron_mouse, h=Baron_human, z=Zeisel

# ---- 构建所有组合 ----
df <- expand.grid(
  method  = methods, 
  dataset = datasets, 
  stringsAsFactors = FALSE
)

# ---- 从环境中安全获取变量 ----
get_times_safe <- function(method, dataset) {
  varname <- paste0(method, "_times", dataset)
  if (exists(varname, envir = .GlobalEnv)) {
    val <- get(varname, envir = .GlobalEnv)
    if (is.numeric(val) && length(val) > 0) return(val)
  }
  return(NA_real_)
}

# ---- 读取时间 ----
df$times <- sapply(seq_len(nrow(df)), function(i) {
  get_times_safe(df$method[i], df$dataset[i])
})

# ---- 创建 fmt（保证不会缺失） ----
df$fmt <- ifelse(
  is.na(df$times),
  "---",
  sprintf("%.3f", df$times)
)

# ---- 转宽表 ----
table_df <- df %>%
  dplyr::select(method, dataset, fmt) %>%
  tidyr::pivot_wider(names_from = dataset, values_from = fmt)

# ---- 美化方法名 ----
table_df$method <- stringr::str_to_title(table_df$method)

# ---- 改列名 ----
table_df <- table_df %>%
  dplyr::rename(
    Method          = method,
    `Baron (Mouse)` = m,
    `Baron (Human)` = h,
    Zeisel          = z
  )

# ---- 生成 LaTeX ----
latex <- paste(
  "\\begin{table}[H]",
  "\\centering",
  "\\caption{Runtime comparison of annotation methods across three datasets.}",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Method & Baron (Mouse) & Baron (Human) & Zeisel \\\\",
  "\\midrule",
  sep = "\n"
)

for (i in seq_len(nrow(table_df))) {
  row <- table_df[i, ]
  latex_row <- sprintf(
    "%s & %s & %s & %s \\\\",
    row$Method,
    row$`Baron (Mouse)`,
    row$`Baron (Human)`,
    row$Zeisel
  )
  latex <- paste(latex, latex_row, sep = "\n")
}

latex <- paste(
  latex,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}",
  sep = "\n"
)

cat(latex)




