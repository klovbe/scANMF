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
# Load dataset human
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

# Load dataset mouse
Romanov <- RomanovBrainData()
unique(colData(Romanov)$'level1 class')
Romanov_label_filtered <- colData(Romanov)$'level1 class'
Romanov <- genem2h(Romanov)

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

Darmanis_label_filtered <- gsub("oligodendrocytes", "Oligodendrocytes", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("OPC", "Oligodendrocyte precursor cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("microglia", "Microglial cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("endothelial", "Endothelial cells", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("neurons", "Neurons", Darmanis_label_filtered)
Darmanis_label_filtered <- gsub("astrocytes", "Astrocytes", Darmanis_label_filtered)

Romanov_label_filtered <- gsub("ependymal", "Ependymal", Romanov_label_filtered)
Romanov_label_filtered <- gsub("vsm", "Smooth Muscle", Romanov_label_filtered)
Romanov_label_filtered <- gsub("neurons", "Neurons", Romanov_label_filtered)
Romanov_label_filtered <- gsub("oligos", "Oligodendrocytes", Romanov_label_filtered)
Romanov_label_filtered <- gsub("microglia", "Microglial cells", Romanov_label_filtered)
Romanov_label_filtered <- gsub("endothelial", "Endothelial cells", Romanov_label_filtered)
Romanov_label_filtered <- gsub("astrocytes", "Astrocytes", Romanov_label_filtered)


checkGeneSymbols <- function(x, ...) {
  data.frame(x, Suggested.Symbol = x, Approved = TRUE, stringsAsFactors = FALSE)
}

source("/Users/weilaichi/Downloads/R/sc-type-master/R/gene_sets_prepare.R")
source("/Users/weilaichi/Downloads/R/sc-type-master/R/sctype_score_.R")
# DB file
library(openxlsx)
tissue = "Brain" # e.g. Immune system,brain,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
db_ <- "/Users/weilaichi/Downloads/R/sc-type-master/ScTypeDB_full.xlsx";
gs_list = gene_sets_prepare(db_, tissue)
gs = gs_list$gs_positive
gs$`Cancer cells` <- NULL
gs$`Cancer stem cells` <- NULL
gs_list1 = gene_sets_prepare(db_, "Hippocampus")
gs1 = gs_list1$gs_positive
gs_new <- c(gs, gs1[c('Ependymal',"Smooth Muscle")])


#------------------------------------------------------------------
#darmanis
#------------------------------
z <- unique(intersect(union(toupper(highly_variable_genes_Romanov),toupper(highly_variable_genes_darmanis)),intersect(toupper(rownames(expression_matrix_Romanov)), toupper(rownames(expression_matrix_Darmanis)))))
# assign cell types
matrix_from_list <- do.call(cbind, gs_new)
celltypes <- colnames(matrix_from_list)
vector<- c(matrix_from_list)
unique_vector <- unique(vector)
unique_vector <- unique(intersect(unique_vector,intersect(toupper(rownames(expression_matrix_Romanov)), toupper(rownames(expression_matrix_Darmanis)))))
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

# 高变基因标记基因
#基因表达矩阵与高变基因标记基因在高边基因内部的部分
rownames(expression_matrix_Darmanis) <- toupper(rownames(expression_matrix_Darmanis))
rownames(expression_matrix_Romanov) <- toupper(rownames(expression_matrix_Romanov))
expression_matrix_common_Darmanis <- expression_matrix_Darmanis[union_vec, , drop = FALSE]
expression_matrix_common_Romanov <- expression_matrix_Romanov[union_vec, , drop = FALSE]
expression_matrix_common <- cbind(expression_matrix_common_Darmanis, expression_matrix_common_Romanov)
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
num_cells_Darmanis <- length(Darmanis_label_filtered)
num_cells_Romanov <- length(Romanov_label_filtered)
num_cells <- num_cells_Romanovnum_cells_Darmanis
#——————————————————————————————————
#####################
#以Romanov预测darmanis
############################
#-----------------------------
#scANMF
#——————————————————————————————————
#标签信息
unique_labels <- unique(Romanov_label_filtered)
num_labels <- length(unique_labels)
expression_matrix_common <- cbind(expression_matrix_common_Romanov, expression_matrix_common_Darmanis)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:length(Romanov_label_filtered)) {
  labels <- Romanov_label_filtered[i]
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
#标记基因对应细胞类型标签
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
#KNNGauss
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
  Q <- ((Q  t(Q)) / 2)
}

Q <- (Q  t(Q)) / 2   # 再次确保对称
D <- diag(rowSums(Q))
L <- D - Q             # 未归一化拉普拉斯
# 若需对称归一化，可用：
# D_inv_sqrt <- diag(1 / sqrt(diag(D)  1e-10))
# L <- diag(n) - D_inv_sqrt %*% Q %*% D_inv_sqrt



# alpha<-1000
# beta<-20
# lambda<-100

alpha<-1000
beta<-30
lambda<-300

result<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
cell_type_matrix  <-result[[2]]
marker<-result[[1]]


library(scales)
marker_brain <- marker
Markers_brain <- Markers
Markers_brain <- Markers_brain[rowSums(Markers_brain != 0) > 0, ]
# 保留不全为0的行
Markers_brain <- 1 - Markers_brain

# 只保留 gene_Muraro 中和 Markers_brain 具有相同行名的行
facotr_brain <- marker_brain[rownames(marker_brain) %in% rownames(Markers_brain), , drop = FALSE]
colnames(facotr_brain)<-colnames(labelMatrix)
facotr_brain <- facotr_brain[rownames(Markers_brain),]

# 输入：
# W:   genes x K (matrix)
# Markers_brain: genes x C (binary matrix; 1 = marker of that cell type)
# gene order应当对齐：rownames(W) == rownames(Markers_brain)
marker_threshold <- 5   # 你可以调成 2、5 等

# 每个细胞类型的 marker 数
marker_count <- colSums(Markers_brain == 1)

# 保留满足条件的 cell types
keep_ct <- names(marker_count[marker_count >= marker_threshold])

cat("保留的细胞类型：\n")
print(keep_ct)

cat("被删除的细胞类型：\n")
print(setdiff(colnames(Markers_brain), keep_ct))

# 筛掉 marker 太少的 cell types
Marker <- Markers_brain[, keep_ct, drop = FALSE]
W <- facotr_brain[, keep_ct, drop = FALSE]

save(Marker,W, file = 'pheat_rom2darm.RData')
U<-W

library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(viridis)

##########################################################
# 1. 修正后的函数：只计算，不画图
##########################################################

evaluate_marker_metrics_weighted <- function(U, Marker){
  
  k_star <- apply(U, 1, which.max)
  marker_sets <- apply(Marker, 1, function(x) which(x == 1))
  
  K <- ncol(Marker)
  p <- ncol(U)
  
  # weight matrix
  W_marker <- Matrix(0, nrow = nrow(Marker), ncol = K,
                     dimnames = dimnames(Marker))
  for(i in 1:nrow(Marker)){
    Cs <- marker_sets[[i]]
    if(length(Cs) > 0){
      W_marker[i, Cs] <- 1 / length(Cs)
    }
  }
  
  # weighted counts n_{c→k}
  n_c2k <- matrix(0, nrow = K, ncol = p,
                  dimnames = list(colnames(Marker), colnames(U)))
  for(c in 1:K){
    for(k in 1:p){
      n_c2k[c,k] <- sum(W_marker[k_star == k, c])
    }
  }
  
  denom_c <- colSums(W_marker)
  P_c2k <- sweep(n_c2k, 1, denom_c, "/")
  
  factor_map <- apply(n_c2k, 2, which.max)
  
  correct_vec <- rep(FALSE, nrow(U))
  for(i in 1:nrow(U)){
    Cs <- marker_sets[[i]]
    if(length(Cs) > 0){
      correct_vec[i] <- factor_map[k_star[i]] %in% Cs
    }
  }
  
  marker_accuracy <- mean(correct_vec[rowSums(Marker)>0])
  
  return(list(
    marker_accuracy = marker_accuracy,
    P_c2k = P_c2k,
    factor_map = factor_map,
    n_c2k = n_c2k
  ))
}

##########################################################
# 2. 调用函数，拿到结果
##########################################################

res <- evaluate_marker_metrics_weighted(U, Marker)
P_c2k <- res$P_c2k
factor_map <- res$factor_map

cat("Marker Accuracy =", res$marker_accuracy, "\n")

##########################################################
# 3. 现在正式画图（不放在函数里）
##########################################################

mat <- P_c2k


library(RColorBrewer)

# color function
col_fun <- colorRamp2(
    breaks = seq(0, max(mat), length.out = 5),
    colors = viridis(5)
    )

# 2. 为 Assigned cell types 定义离散颜色
cell_types <- unique(colnames(Marker))
assign_colors <- c(
  "Immune system cells"             = "#D55E00",  # orange/red
  "Myelinating Schwann cells"       = "#009E73",  # bluish green
  "Neural Progenitor cells"         = "#56B4E9",  # sky blue
  "Neural stem cells"               = "#0072B2",  # blue
  "Neuroblasts"                     = "#4DBBD5",  # teal/cyan
  "Neuroepithelial cells"           = "#00A087",  # dark teal
  "Non myelinating Schwann cells"   = "#008280",  # greenish cyan
  "Oligodendrocyte precursor cells" = "#7E6148",  # brown
  "Radial glial cells"              = "#CC79A7",  # magenta
  "Schwann precursor cells"         = "#009988",  # turquoise
  "Tanycytes"                       = "#F0E442",  # yellow
  "Oligodendrocytes"                = "#6A51A3",  # deep purple
  "Astrocytes"                      = "#7570B3",  # bluish purple
  "Ependymal"                       = "#A6761D",  # deep brown
  "Microglial cells"                = "#E64B35",  # red-orange
  "Smooth Muscle"                   = "#B2182B",  # deep red
  "Endothelial cells"               = "#EF8A62",  # salmon
  "Neurons"                         = "#1B9E77"   # green
)
cell_types <- colnames(Markers)

# 确保完全匹配
assign_colors <- assign_colors[cell_types]

# 3. Top annotation
col_ha <- HeatmapAnnotation(
  'Assigned Cell Type' = factor(colnames(Marker)[factor_map], levels = cell_types),
  col = list(Assigned = assign_colors),
  annotation_height = unit(4, "mm")
)


h <- Heatmap(
     mat,
     name = "Proportion",
     col = col_fun,
     cluster_rows = FALSE,
     cluster_columns = FALSE,
     row_names_side = "left",
     column_names_side = "bottom",
     column_title = "Factors",
     row_title = "Cell Types",
     heatmap_legend_param = list(
         title = "Marker proportion"
       ),
     bottom_annotation = col_ha
   )
# 输出到 PDF（矢量图）



pdf("heatmap.pdf", width = 12, height = 8)

grid::grid.newpage()     # 新建页
draw(h, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     newpage = FALSE,    # 很重要！否则会再开启一个空页面
     merge_legend = TRUE,
     padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()



