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
true_labels <- c(mouse_label_filtered, human_label_filtered)
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
# 计算注释准确率
true_labels<-human_label_filtered
annotation_labels<-annotation_labels[(num_cells_mouse+1):length(annotation_labels)]
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
cross_scanmf_m2h <- metric_cal(true_labels, annotation_labels)

######################
#SingleR 0.8982456
library(SingleR)
pred<- SingleR(test=expression_matrix_human, ref=expression_matrix_mouse, labels=mouse_label_filtered)
pred_singler <- pred$labels
accuracy <- sum(human_label_filtered == pred_singler) / length(human_label_filtered)
# 输出准确率
print(accuracy)
table(pred_singler, human_label_filtered)

cross_singler_m2h <- metric_cal(human_label_filtered, pred_singler)

#####scpred 0.3929825
# Load dataset mouse
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
acc <- mean(pred_scpred == query_label)
print(acc)
table(pred_scpred, query_label)

cross_scpred_m2h <- metric_cal(query_label, pred_scpred)


#######################################
#以human预测mouse
#######################################

#-----------------------------
#scANMF
#——————————————————————————————————
#标签信息
retain_ratio <- 0.2  # 例如只保留 50% 细胞的标签
set.seed(123)
n <- length(human_label_filtered)
retain_mask <- rep(FALSE, n)
retain_mask[sample(seq_len(n), retain_num)] <- TRUE
# retain_mask[which(human_label_filtered %in% c("Mast cells", "Epsilon cells", "Peri-islet Schwann cells"))]<-FALSE
human_label_filtered_part <- human_label_filtered[retain_mask]
expression_matrix_common_human_part <- expression_matrix_common_human[,retain_mask]
unique_labels <- unique(human_label_filtered_part)
num_labels <- length(unique_labels)
expression_matrix_common <- cbind(expression_matrix_common_human_part, expression_matrix_common_mouse)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
label_matrix <- matrix(0, nrow = ncol(expression_matrix_common), ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:ncol(expression_matrix_common_human_part)) {
  labels <- human_label_filtered[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}


colnames(label_matrix)<-unique_labels
#标签和标记基因对应的细胞类型的并集
unique_labellist <- unique(c(colnames(markers_common),colnames(label_matrix)))
#标签以外的细胞类型
a<- setdiff(unique_labellist, colnames(label_matrix))
#额外的细胞类型矩阵
label_matrix1 <- matrix(0, nrow = ncol(expression_matrix_common), ncol =length(a) )
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


alpha<-1000000
beta<-1
lambda<-20



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
true_labels <- c(human_label_filtered_part, mouse_label_filtered)
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
# 计算注释准确率
annotation_labels<-annotation_labels[(length(human_label_filtered_part)+1):length(annotation_labels)]
true_labels <- mouse_label_filtered
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
cross_scanmf_h2m <- metric_cal(mouse_label_filtered, annotation_labels)




retain_ratio <- 0.2  # 例如只保留 50% 细胞的标签
set.seed(123)
n <- length(human_label_filtered)
retain_num <- round(n * retain_ratio)
retain_mask <- rep(FALSE, n)
retain_mask[sample(seq_len(n), retain_num)] <- TRUE
# retain_mask[which(human_label_filtered %in% c("Mast cells", "Epsilon cells", "Peri-islet Schwann cells"))]<-TRUE
human_label_filtered_part <- human_label_filtered[retain_mask]
expression_matrix_common_human_part <- expression_matrix_common_human[,retain_mask]
unique_labels <- unique(human_label_filtered_part)
num_labels <- length(unique_labels)
expression_matrix_common <- cbind(expression_matrix_common_human_part, expression_matrix_common_mouse)


#singleR
expression_matrix_human_part <- expression_matrix_human[,retain_mask]
pred<- SingleR(test=expression_matrix_mouse, ref=expression_matrix_human_part, labels=human_label_filtered_part)
pred_singler <- pred$labels
true_labels <- mouse_label_filtered
accuracy <- sum(true_labels == pred_singler) / length(true_labels)
# 输出准确率
print(accuracy)
table(pred_singler, mouse_label_filtered)
cross_singler_h2m <- metric_cal(mouse_label_filtered, pred_singler)

#####scpred 0.3929825
# Load dataset mouse
retain_ratio <- 0.2  # 例如只保留 50% 细胞的标签
set.seed(123)
n <- length(human_label_filtered)
retain_num <- round(n * retain_ratio)
retain_mask <- rep(FALSE, n)
retain_mask[sample(seq_len(n), retain_num)] <- TRUE
retain_mask[which(human_label_filtered %in% c("Mast cells", "Epsilon cells", "Peri-islet Schwann cells"))]<-TRUE
human_filtered_part <- human_filtered[,retain_mask]
ref <- GetAssayData(human_filtered_part, assay = "RNA",layer = "counts")
query<-GetAssayData(mouse_filtered, assay = "RNA",layer = "counts")
human_label_filtered_part <- human_label_filtered[retain_mask]
ref_label <- human_label_filtered_part
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
acc <- mean(pred_scpred == query_label)
print(acc)
table(pred_scpred, query_label)

cross_scpred_h2m <- metric_cal(query_label, pred_scpred)

#---------------
##sccatch
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

## 检查结果
table(pred_aligned, truth_aligned)
acc <- mean(pred_aligned == truth_aligned)

# 混淆矩阵（可选）
conf <- table(pred = pred_aligned, truth = truth_aligned)
cross_sccatch_h <- metric_cal(human_label_filtered, pred_aligned)


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

## 检查结果
table(pred_aligned, truth_aligned)
acc <- mean(pred_aligned == truth_aligned)

# 混淆矩阵（可选）
conf <- table(pred = pred_aligned, truth = truth_aligned)
cross_sccatch_m <- metric_cal(mouse_label_filtered, pred_aligned)

#sctype
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

## 检查结果
conf <- table(pred_aligned, truth_aligned)
acc <- mean(pred_aligned == truth_aligned)

cross_sctype_m <- metric_cal(mouse_label_filtered, pred_aligned)



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

## 检查结果
conf <- table(pred_aligned, truth_aligned)
acc <- mean(pred_aligned == truth_aligned)

cross_sctype_h <- metric_cal(human_label_filtered, pred_aligned)
cross_scpred_h2m <- list()

save(cross_scanmf_h2m, cross_scanmf_m2h, cross_scpred_h2m, cross_scpred_m2h, cross_singler_h2m, cross_singler_m2h, cross_sccatch_h, cross_sccatch_m, cross_sctype_h, cross_sctype_m, file = "baron_cross.RData")

load("baron_cross.RData")



