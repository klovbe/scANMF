# ===== 依赖 =====
suppressPackageStartupMessages({
  library(data.table)
  library(SeuratObject)
  library(stringr)
  library(Matrix)
  library(Seurat)
  library(scCATCH)
  library(SeuratObject)
  library(dplyr)
  library(FNN)
  library(scPred)
  library(SingleCellExperiment)
  library(scRNAseq)
  library(ggplot2)
  library(scales)
  library(SingleR)
  library(openxlsx)
})

run_pancreas <- function(result_sc,
                                 method = "scCATCH",
                                 seurat_min_cells = 3,
                                 seurat_min_features = 200,
                                 norm_method = "LogNormalize",
                                 n_top_var_features = 2000,
                                 resolution = 0.5,
                                 only_pos_markers = TRUE,
                                 min_pct = 0.25,
                                 logfc_threshold = 0.25,
                                 seed = 1234) {
  # stopifnot(inherits(result_sc$te, "dgCMatrix"))
  # stopifnot(length(result_sc$label) == ncol(result_sc$te))
  # set.seed(seed)

  
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
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE, nfeatures = n_top_var_features)
  seu <- Seurat::ScaleData(seu, features = rownames(seu), verbose = FALSE)
  seu <- Seurat::RunPCA(seu, features = VariableFeatures(seu), verbose = FALSE)
  seu <- Seurat::FindNeighbors(seu, verbose = FALSE)
  seu <- Seurat::FindClusters(seu, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:10, graph = NULL, nn.name = NULL, features = NULL)
  DimPlot(seu, reduction = "umap", label = TRUE)
  expression_matrix_pan <- GetAssayData(seu, assay = "RNA")
  
  if (method == 'singleR'){
    #singleR
    x <- expression_matrix_pan[,1:10600]
    y <- expression_matrix_pan[,10601:length(truth)]
    true_labels <- truth[10601:length(truth)]
    pred <- SingleR(y,x,labels = truth[1:10600])
    accuracy <- sum(true_labels == pred$labels) / length(true_labels)
    # 输出准确率
    print(accuracy)
    table(pred$labels, true_labels)  
    result_singleR_forward <-   metric_cal(true_labels, pred$labels)
    
    y <- expression_matrix_pan[,1:10600]
    x <- expression_matrix_pan[,10601:length(truth)]
    true_labels <- truth[1:10600]
    pred <- SingleR(y,x,labels = truth[10601:length(truth)])
    accuracy <- sum(true_labels == pred$labels) / length(true_labels)
    # 输出准确率
    print(accuracy)
    table(pred$labels, true_labels)  
    result_singleR_backward <-   metric_cal(true_labels, pred$labels)    
    return(list(result_singleR_forward, result_singleR_backward))
  }  

  else if (method == 'scpred'){
    #scPred
    query_label <- truth[10601:length(truth)]
    ref_label <- truth[1:10600]
    
    ref <- CreateSeuratObject(
      counts = expr[,1:10600],
      project = "pancreas"
    )
    DefaultAssay(ref) <- "RNA"
    ref <- SetAssayData(ref, assay = "RNA", layer = "data", new.data = expr)
    ref <- FindVariableFeatures(ref, verbose = FALSE)
    ref <- ScaleData(ref, features = rownames(ref), verbose = FALSE)
    ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 10, reduction.key = "pca_", verbose = FALSE, approx = FALSE)
    ref@meta.data$cell_type <- truth[match(colnames(ref), cells)]
    DefaultAssay(ref) <- "RNA"
    
    
    query <- CreateSeuratObject(counts = expr[,10601:length(truth)])
    query <- NormalizeData(query)
    DefaultAssay(query) <- "RNA"
    ref <- getFeatureSpace(ref, "cell_type")
    ref <- trainModel(ref)
    query <- scPredict(query, ref)
    # DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
    pred_scpred <- query@meta.data$scpred_prediction  
    result_scpred_forward <- metric_cal(query_label, pred_scpred)
    table(query_label, pred_scpred)
    
    
    
    ref_label <- truth[10601:length(truth)]
    query_label <- truth[1:10600]
    
    ref <- CreateSeuratObject(
      counts = expr[,10601:length(truth)],
      project = "pancreas"
    )
    DefaultAssay(ref) <- "RNA"
    ref <- SetAssayData(ref, assay = "RNA", layer = "data", new.data = expr)
    ref <- FindVariableFeatures(ref, verbose = FALSE)
    ref <- ScaleData(ref, features = rownames(ref), verbose = FALSE)
    ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = 10, reduction.key = "pca_", verbose = FALSE, approx = FALSE)
    ref@meta.data$cell_type <- truth[match(colnames(ref), cells)]
    DefaultAssay(ref) <- "RNA"
    
    
    query <- CreateSeuratObject(counts = expr[,1:10600])
    query <- NormalizeData(query)
    DefaultAssay(query) <- "RNA"
    ref <- getFeatureSpace(ref, "cell_type")
    ref <- trainModel(ref)
    query <- scPredict(query, ref)
    # DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
    pred_scpred <- query@meta.data$scpred_prediction  
    
    
    result_scpred_backward <- metric_cal(query_label, pred_scpred)
    table(query_label, pred_scpred)
    return(list(result_scpred_forward, result_scpred_backward))
  }

  else if (method == 'scanmf'){
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
    
    max_col_indices <- apply(cell_type_matrix, 1, which.max)
    ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
    normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix), ncol = ncol(cell_type_matrix))
    for (j in seq_along(max_col_indices)) {
      normalized_matrix[j, max_col_indices[j]] <- 1
    }
    cell_type_matrix <- normalized_matrix
    colnames(cell_type_matrix)<-colnames(labelMatrix)
    annotation <-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
    # 计算注释准确率
    true_labels<-truth_mapped[!retain_mask]
    annotation_labels<-annotation[!retain_mask]
    accuracy <- sum(true_labels == annotation_labels) / length(true_labels)
    # 输出准确率
    print(accuracy)
    table(annotation_labels, true_labels)
    result_scanmf_forward <- metric_cal(true_labels, annotation_labels)
    
    
    #backward
    retain_mask <- rep(FALSE, length(truth_mapped))
    retain_mask[10601:num_cells] <- TRUE
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
    
    max_col_indices <- apply(cell_type_matrix, 1, which.max)
    ##创建一个新的矩阵，只保留每一列最大值所在位置的值，并将这些值标准化为1
    normalized_matrix <- matrix(0, nrow = nrow(cell_type_matrix), ncol = ncol(cell_type_matrix))
    for (j in seq_along(max_col_indices)) {
      normalized_matrix[j, max_col_indices[j]] <- 1
    }
    cell_type_matrix <- normalized_matrix
    colnames(cell_type_matrix)<-colnames(labelMatrix)
    annotation <-apply(cell_type_matrix, 1, function(row) colnames(cell_type_matrix)[which(row == 1)])
    # 计算注释准确率
    true_labels<-truth_mapped[!retain_mask]
    annotation_labels<-annotation[!retain_mask]
    # 建立映射关系（统一到 markers_common 的命名体系）
    mapping <- c(
      "macrophage" = "Immune system cells"
    )
    true_labels1 <- ifelse(true_labels %in% names(mapping), mapping[true_labels], true_labels)
    accuracy <- sum(true_labels1 == annotation_labels) / length(true_labels1)
    # 输出准确率
    print(accuracy)
    table(annotation_labels, true_labels1)
    result_scanmf_backward <- metric_cal(true_labels1, annotation_labels)
    return(list(result_scanmf_forward, result_scanmf_backward))
  }

  else if (method == "scCATCH"){
    obj <- scCATCH::createscCATCH(
      data    = as.matrix(GetAssayData(seu, assay="RNA", layer="data")),
      cluster = as.character(Idents(seu))
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
    seu$scCATCH_celltype <- unname(cl_map[as.character(Idents(seu))])
    pred_celltype <- seu$scCATCH_celltype
    
    cat("per-cluster 注释条目：", nrow(ct), "\n",
        "逐细胞 NA 比例：", mean(is.na(seu$scCATCH_celltype)), "\n")
    
    
    truth_aligned <- seu$Celltype2_truth
    map_unique <- data.frame(
      pred = c(
        "Acinar Cell",
        "Alpha Cell",
        "Beta Cell",
        "Beta Cell, Stem Cell",
        "Delta Cell",
        "PP Cell",
        "Mesenchymal Cell, Stem Cell"
      ),
      mapped_truth = c(
        "acinar",
        "alpha",
        "beta",
        "beta",
        "delta",
        "PP",
        "mesenchymal"
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
    
    ## 检查结果
    table(pred_aligned, truth_aligned)
    acc <- mean(pred_aligned == truth_aligned)
    
    # 混淆矩阵（可选）
    conf <- table(pred = pred_aligned, truth = truth_aligned)
    
    list(
      seurat_obj     = seu,
      scCATCH_obj    = obj,
      pred  = pred_aligned,         # 每个细胞的预测标签
      truth = truth_aligned,         # 每个细胞的真值
      accuracy       = acc,
      confusion      = conf
    )    
  }  
  else if (method == "ScType"){
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
    
    min_score <- NULL
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
    truth_aligned <- seu$Celltype2_truth
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
    truth_aligned <- ifelse(truth_aligned %in% names(mapping), mapping[truth_aligned], truth_aligned)
    pred_aligned <- labels_cluster
    pred_aligned[is.na(pred_aligned)] <- "unknown"
    ## 去掉名字属性（保持向量干净）
    pred_aligned <- unname(pred_aligned)
    
    ## 检查结果
    conf <- table(pred_aligned, truth_aligned)
    acc <- mean(pred_aligned == truth_aligned)
    
    list(
      seurat_obj     = seu,
      sctype_sc  = es.max,
      pred  = pred_aligned,         # 每个细胞的预测标签
      truth = truth_aligned,         # 每个细胞的真值
      accuracy       = acc,
      confusion      = conf
    ) 
  }

}

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

# 运行scCATCH：
out_catch <- run_pancreas(result_sc, method = "scCATCH")
#case1: 
metric_scCATCH <- metric_cal(out_catch$truth[1:10600],out_catch$pred[1:10600])

#case2: 
metric_scCATCH1 <- metric_cal(out_catch$truth[10601:length(out_catch$truth)],out_catch$pred[10601:length(out_catch$truth)])

# 运行ScType：
out_sctype <- run_pancreas(result_sc, method = "ScType")
#case1: 
metric_sctype <- metric_cal(out_sctype$truth[1:10600],out_sctype$pred[1:10600])

#case2: 
metric_sctype1 <- metric_cal(out_sctype$truth[10601:length(out_sctype$truth)],out_sctype$pred[10601:length(out_sctype$truth)])

#运行scANMF
metric_scanmf <- run_pancreas(result_sc, method = "scanmf")
result_scanmf_forward <- metric_scanmf[[1]]
result_scanmf_backward <- metric_scanmf[[2]]
metric_singleR <- run_pancreas(result_sc, method = "singleR")
metric_scpred <- run_pancreas(result_sc, method = "scpred")

save(result_scanmf_forward, result_scanmf_backward, result_scpred_forward, result_scpred_backward, result_singleR_forward, result_singleR_backward, 
     metric_scCATCH, metric_scCATCH1, metric_sctype, metric_sctype1, file = "cross_pancreas.RData")
