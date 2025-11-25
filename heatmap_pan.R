# 1) 读入 Matrix Market 稀疏矩阵（基因×细胞）
expr <- Matrix::readMM("/Users/weilaichi/Downloads/R/pan_expr.mtx")  # dgCMatrix, genes x cells

# 2) 读入基因和细胞名，并赋值到 dimnames
genes <- fread("/Users/weilaichi/Downloads/R/pan_genes.tsv", header = FALSE)$V1
cells <- fread("/Users/weilaichi/Downloads/R/pan_barcodes.tsv", header = FALSE)$V1
stopifnot(nrow(expr) == length(genes), ncol(expr) == length(cells))
rownames(expr) <- genes
colnames(expr) <- cells

# 3) 读入标签（与列顺序一致）
labels <- fread("/Users/weilaichi/Downloads/R/pan_labels.csv")  # 列: cell, Celltype2, Batch, disease
stopifnot(all(labels$cell == colnames(expr)))
truth_vec <- labels$Celltype2  # 真值标签（字符型）

# ===== 给 scType / scCATCH 使用 =====
# 你的两套流程都能吃 “基因×细胞” 的稀疏矩阵。下面给出一个最小示例：

# (A) 如果你已有 scCATCH 的包装函数 run_scCATCH_with_accuracy()
# 假设它需要一个 result 列表：result$te 为表达矩阵（基因×细胞），result$label 为真值
result_sc <- list(
  te = expr,
  label = truth_vec  # 若你的函数里直接比较字符串标签，OK；若强制 Type1..K，则见注释*
)
expr <- result_sc$te            # genes x cells
truth <- as.character(result_sc$label)
cells <- colnames(expr)
stopifnot(length(truth) == length(cells))
if (!inherits(expr, "CsparseMatrix")) {
  expr <- as(expr, "CsparseMatrix")   # 取代 as(...,"dgCMatrix")
}
# 2) CreateSeuratObject（不 Normalize），把矩阵放入 "data" layer（v5 正确做法）
seu <- Seurat::CreateSeuratObject(
  counts = expr,
  project = "pancreas"
)
DefaultAssay(seu) <- "RNA"
seu <- SetAssayData(seu, assay = "RNA", layer = "data", new.data = expr)
# 这里如果 counts 为空，仍不建议跑 vst；建议改用 dispersion（见上）
# 如一定要 vst，需要 counts：把 raw counts 放入 counts，再 NormalizeData 后 FindVariableFeatures(vst)
# seu <- FindVariableFeatures(seu, assay = "RNA", selection.method = "vst",
#                               nfeatures = 2000, layer = "data")
# 
# 把真值塞进 meta.data，便于后面评估
seu$Celltype2_truth <- truth[match(colnames(seu), cells)]

# --- 标准流程：归一化/找高变/PCA/邻居/聚类 ---
seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE, nfeatures = 2000)
seu <- Seurat::ScaleData(seu, features = rownames(seu), verbose = FALSE)
expression_matrix_pan <- GetAssayData(seu, assay = "RNA")
highly_variable_genes_pan <- VariableFeatures(seu, n=500)
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
union_vec <- toupper(union(highly_variable_genes_pan, unique_vector))
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
expression_matrix_pan<-as.matrix(expression_matrix_pan)

# 高变基因+标记基因
#基因表达矩阵与高变基因+标记基因相交部分
common_rows <- intersect(toupper(rownames(expression_matrix_pan)), toupper(rownames(markers)))
# common_rows 现在包含了两个矩阵中行名相同的行的标识符
common_rows1 <- rownames(expression_matrix_pan)[toupper(rownames(expression_matrix_pan)) %in% common_rows]
expression_matrix_common <- expression_matrix_pan[common_rows1, , drop = FALSE]
rownames(expression_matrix_common)<-toupper(rownames(expression_matrix_common)) 
markers_common<- markers[common_rows, ]
colnames(markers_common)

# 建立映射关系（统一到 markers_common 的命名体系）
mapping <- c(
  "acinar" = "Acinar cells",
  "alpha" = "Alpha cells",
  "beta" = "Beta cells",
  "delta" = "Delta cells",
  "ductal" = "Ductal cells",
  "epsilon" = "Epsilon cells",
  "PP" = "Gamma (PP) cells",
  "endothelial" = "Endothelial cells",
  "mast" = "Mast cells",
  "mesenchymal" = "Mesenchymal cells",
  "PSC" = "Pancreatic stellate cells",
  "schwann" = "Peri-islet Schwann cells",
  "macrophage" ="Immune system cells",
  "t_cell"="Immune system cells"
)

# 替换 truth 中的名称（不在 mapping 的保持原样）
truth_mapped <- ifelse(truth %in% names(mapping), mapping[truth], truth)
# # 查看替换结果
# table(truth, truth_mapped)
#——————————————————————————————————
#标签信息
unique_labels <- unique(truth_mapped)
num_labels <- length(truth_mapped)

# 创建空的标签矩阵，行数为细胞数量，列数为标签类别的数量
num_cells <- length(truth_mapped)
label_matrix <- matrix(0, nrow = num_cells, ncol = length(unique_labels))

# 填充标签矩阵
for (i in 1:num_cells) {
  labels <- truth_mapped[i]
  col_index <- which(unique_labels == labels)
  label_matrix[i, col_index] <- 1
}
# 创建一个示例细胞标签矩阵
n_rows <- num_cells  # 矩阵行数
n_cols <- num_labels   # 矩阵列数

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

#forward
retain_mask <- rep(FALSE, length(truth_mapped))
retain_mask[1:10600] <- TRUE
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

#1000,20,1000:0.9
alpha<-10000
beta<-10
lambda<-10000

result<-scANMF(expression_matrix_common,labelMatrix,ncol(labelMatrix),Markers,alpha,beta,lambda,L,D,Q,50,0.01)
cell_type_matrix  <-result[[2]]
marker <- result[[1]]

library(scales)
marker_pancreas <- marker
Markers_pancreas <- Markers
Markers_pancreas <- Markers_pancreas[rowSums(Markers_pancreas != 0) > 0, ]
# 保留不全为0的行
Markers_pancreas <- 1 - Markers_pancreas

# 只保留 gene_Muraro 中和 Markers_pancreas 具有相同行名的行
factor_pancreas <- marker_pancreas[rownames(marker_pancreas) %in% rownames(Markers_pancreas), , drop = FALSE]
colnames(factor_pancreas)<-colnames(labelMatrix)
factor_pancreas <- factor_pancreas[rownames(Markers_pancreas),]

# 输入：
# W:   genes x K (matrix)
# Markers_pancreas: genes x C (binary matrix; 1 = marker of that cell type)
# gene order应当对齐：rownames(W) == rownames(Markers_pancreas)
marker_threshold <- 5   # 你可以调成 2、5 等

# 每个细胞类型的 marker 数
marker_count <- colSums(Markers_pancreas == 1)

# 保留满足条件的 cell types
keep_ct <- names(marker_count[marker_count >= marker_threshold])

cat("保留的细胞类型：\n")
print(keep_ct)

cat("被删除的细胞类型：\n")
print(setdiff(colnames(Markers_pancreas), keep_ct))

# 筛掉 marker 太少的 cell types
Marker <- Markers_pancreas[, keep_ct, drop = FALSE]
W <- factor_pancreas[, keep_ct, drop = FALSE]
save(Marker,W, file = 'pheat_pancreas.RData')

stopifnot(all(rownames(W) == rownames(Markers_pancreas)))
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



pdf("heatmap_pan.pdf", width = 12, height = 8)

grid::grid.newpage()     # 新建页
draw(h, 
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     newpage = FALSE,    # 很重要！否则会再开启一个空页面
     merge_legend = TRUE,
     padding = unit(c(5, 5, 5, 5), "mm")
)

dev.off()




