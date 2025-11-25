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
Zeisel<- ZeiselBrainData()
unique(colData(Zeisel)$level1class)
Zeisel_label <- colData(Zeisel)$level1class
Zeisel_label_filtered <- Zeisel_label
Zeisel_filtered <- Zeisel_filtered <- genem2h(Zeisel)


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

# Load dSingleR# Load dataset human
Darmanis<- DarmanisBrainData()
Darmanis_label <- colData(Darmanis)$cell.type
unique(Darmanis_label)
# 找到 fetal_quiescent 和 hybrid 的索引
indices_to_remove <- which(Darmanis_label %in% c("fetal_quiescent", "hybrid","fetal_replicating"))
# indices_to_remove <- which(Darmanis_label == "hybrid")
# 删除对应的标签信息
Darmanis_label_filtered <- Darmanis_label[-indices_to_remove]
Darmanis_filtered <- Darmanis[,-indices_to_remove]

Darmanis_filtered <- CreateSeuratObject(counts(Darmanis_filtered),project = "Darmanis",min.cells = 3, min.features = 200) #创建seurat对象
#QC指标
Darmanis_filtered[["percent.mt"]] <- PercentageFeatureSet(Darmanis_filtered, pattern = "^MT-")
Darmanis_filtered <- subset(Darmanis_filtered, subset = percent.mt < 20)
Darmanis_filtered <- NormalizeData(Darmanis_filtered)
Darmanis_filtered <- FindVariableFeatures(Darmanis_filtered)
Darmanis_filtered<- ScaleData(Darmanis_filtered, features = rownames(Darmanis_filtered))
highly_variable_genes_darmanis <- VariableFeatures(Darmanis_filtered, n=1000)
expression_matrix_Darmanis <- GetAssayData(Darmanis_filtered, assay = "RNA")

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

Darmanis_label_filtered <- gsub("oligodendrocytes", "Oligodendrocytes", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("OPC", "Oligodendrocyte precursor cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("microglia", "Microglial cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("endothelial", "Endothelial cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("neurons", "Neurons", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("astrocytes", "Astrocytes", Darmanis_label_filtered)



z <- unique(intersect(union(toupper(highly_variable_genes_darmanis),toupper(highly_variable_genes_zeisel)),intersect(toupper(rownames(expression_matrix_Darmanis)), toupper(rownames(expression_matrix_Zeisel)))))
load("Brain_marker.RData")
gs$`Cancer cells` <- NULL
gs$`Cancer stem cells` <- NULL
# assign cell types
matrix_from_list <- do.call(cbind, gs)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
unique_vector <- unique(intersect(unique_vector,intersect(toupper(rownames(expression_matrix_Darmanis)), toupper(rownames(expression_matrix_Zeisel)))))
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
rownames(expression_matrix_Zeisel) <- toupper(rownames(expression_matrix_Zeisel))
rownames(expression_matrix_Darmanis) <- toupper(rownames(expression_matrix_Darmanis))
expression_matrix_common_zeisel <- expression_matrix_Zeisel[union_vec, , drop = FALSE]
expression_matrix_common_darmanis <- expression_matrix_Darmanis[union_vec, , drop = FALSE]
expression_matrix_common <- cbind(expression_matrix_common_zeisel, expression_matrix_common_darmanis)
markers_common<- markers
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


#——————————————————————————————————
#标签信息
unique_labels <- unique(Zeisel_label_filtered)
num_labels <- length(unique_labels)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells_zeisel <- length(Zeisel_label_filtered)
num_cells_darmanis <- length(Darmanis_label_filtered)
num_cells <- num_cells_darmanis+num_cells_zeisel
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells_zeisel) {
  labels <- Zeisel_label_filtered[i]
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
alpha<-2000
beta<-50
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
# 计算注释准确率
true_labels<-Darmanis_label_filtered
annotation_labels<-annotation_labels[(num_cells_zeisel+1):length(annotation_labels)]
annotation_labels <- gsub("Glutamatergic neurons", "Neurons",annotation_labels)
annotation_labels <- gsub("GABAergic neurons","Neurons" ,annotation_labels)
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
cross_scanmf_z2d <- metric_cal(true_labels, annotation_labels)

######################
#SingleR 0.8982456
library(SingleR)
pred<- SingleR(test=expression_matrix_Darmanis, ref=expression_matrix_Zeisel, labels=Zeisel_label_filtered)
pred_singler <- pred$labels
pred_singler <- gsub("Glutamatergic neurons", "Neurons",pred_singler)
pred_singler <- gsub("GABAergic neurons","Neurons" ,pred_singler)
accuracy <- sum(Darmanis_label_filtered == pred_singler) / length(Darmanis_label_filtered)
# 输出准确率
print(accuracy)
table(pred_singler, Darmanis_label_filtered)

cross_singler_z2d <- metric_cal(Darmanis_label_filtered, pred_singler)

#####scpred 0.3929825
# Load dataset mouse
indices_to_remove <- which(Darmanis_label %in% c("fetal_quiescent", "hybrid","fetal_replicating"))
Darmanis_filtered <- Darmanis[,-indices_to_remove]
Zeisel_filtered <- Zeisel
ref <- counts(Zeisel_filtered)
query<-counts(Darmanis_filtered)
ref_label <- Zeisel_label_filtered
query_label <- Darmanis_label_filtered

ref <- CreateSeuratObject(ref, project = "ref", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = percent.mt < 20)
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

cross_scpred_z2d <- metric_cal(query_label, pred_scpred)

#####################
#以Darmanis预测Zeisel
############################
library(SingleR)
pred<- SingleR(test=expression_matrix_Zeisel, ref=expression_matrix_Darmanis, labels=Darmanis_label_filtered)
pred_singler <- pred$labels
true_labels <- Zeisel_label_filtered
true_labels <- gsub("Glutamatergic neurons", "Neurons",true_labels)
true_labels <- gsub("GABAergic neurons","Neurons" ,true_labels)
accuracy <- sum(true_labels == pred_singler) / length(true_labels)
# 输出准确率
print(accuracy)
table(pred_singler, Zeisel_label_filtered)
cross_singler_d2z <- metric_cal(Zeisel_label_filtered, pred_singler)

#####scpred
# Load dataset mouse
indices_to_remove <- which(Darmanis_label %in% c("fetal_quiescent", "hybrid","fetal_replicating"))
Darmanis_filtered <- Darmanis[,-indices_to_remove]
Zeisel_filtered <- Zeisel
ref <- counts(Darmanis_filtered)
query<-counts(Zeisel_filtered)
ref_label <- Darmanis_label_filtered
query_label <- Zeisel_label_filtered

ref <- CreateSeuratObject(ref, project = "ref", min.cells = 3, min.features = 200)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
ref <- subset(ref, subset = percent.mt < 20)
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
acc <- mean(pred_scpred == true_labels)
print(acc)
table(pred_scpred, true_labels)
cross_scpred_d2z <- metric_cal(Zeisel_label_filtered, pred_scpred)
#-----------------------------
#scANMF
#——————————————————————————————————
#标签信息
unique_labels <- unique(Darmanis_label_filtered)
num_labels <- length(unique_labels)
expression_matrix_common <- cbind(expression_matrix_common_darmanis, expression_matrix_common_zeisel)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells_darmanis) {
  labels <- Darmanis_label_filtered[i]
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



alpha<-2000
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
# 计算注释准确率
annotation_labels<-annotation_labels[(num_cells_darmanis+1):length(annotation_labels)]
accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
# 输出准确率
print(accuracy)
table(annotation_labels, true_labels)
cross_scanmf_d2z <- metric_cal(Zeisel_label_filtered, annotation_labels)

save(cross_scanmf_d2z, cross_scanmf_z2d, cross_scpred_d2z, cross_scpred_z2d, cross_singler_d2z, cross_singler_z2d, file = "darmanis_zeisel.RData")
