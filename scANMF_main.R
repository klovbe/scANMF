########最原始的：
#基于单细胞转录组测序数据的细胞注释算法研究主代码
#X-基因表达矩阵，A标签矩阵，Marker标记基因矩阵，alpha-标记权重，beta-几何权重，lambda-标签权重，iteration-迭代次数，sigma-阈值
scANMF<-function(X,A,rank,Marker,alpha,beta,lambda,L,D,Q,iteration,sigma,init='runif'){
  ##initial W & H
  #kkk = initial_NMF(X,rank,"nndsvd")（初始方法可选）
  set.seed(123)
  num_rows_W <- nrow(X)
  num_cols_W <- rank
  num_rows_H <- ncol(X)
  num_cols_H <- rank
  # 生成随机初始化的矩阵 U 和 Z
  if(init == 'runif'){
    W <- matrix(runif(num_rows_W * num_cols_W), nrow = num_rows_W, ncol = num_cols_W)
    H <- matrix(runif(num_rows_H * num_cols_H), nrow = num_rows_H, ncol = num_cols_H)
  } 

  else{
    init <- nndsvd_init(X, rank, variant = "nndsvda")  # 或 "nndsvd"/"nndsvdar"
    W <- init$W
    H <- init$H  
  }

  # library(NMF)
  # init <- nmf(X, rank, method = "nndsvd", .options = "t")
  # W <- basis(init)
  # H <- coef(init)
  
  
  MM <- alpha*Marker
  zero_rows <- rowSums(MM == 0) == ncol(MM)
  MM[zero_rows, ] <- 50
  
  AA <- lambda*A
  # zero_rows <- rowSums(AA == 0) == ncol(AA)
  # AA[zero_rows, ] <- lambda/10  
  
  
  # 打印生成的矩阵
  
  #计算目标函数
  # objective<-norm(X-W%*%t(H),"F")^2+alpha * sum(W * Marker)+beta*sum(diag(t(H) %*% L %*%H))+lambda*sum(H*A)
  cha<-1
  Ti <-0
  objective<-0
  tt<-numeric(length = iteration)
  while (abs(cha)>sigma && Ti<=iteration){
    # print(objective)
    # print(norm(X-W%*%t(H),"F")^2)
    # print(alpha * sum(W * Marker))
    # print(beta*sum(diag(t(H) %*% L %*%H)))
    # print(lambda*sum(H*A))
    # 计算 X V
    XH <- X %*% H
    
    # 计算 UV^T V
    WHtH <- W %*% t(H) %*% H
    # 计算更新后的 U
    W= W * (2 * XH) / ((2 * WHtH + MM))  # 添加一个小的常数避免除零错误
    # 计算更新后的 V
    H=H * (2*(t(X) %*% W  + 2*beta *  Q %*% H)) / (2*(H %*%t(W) %*% W +2* beta *  D %*% H)+AA)
    #H = H * (t(X) %*% W + beta * Q %*% H) / ( H %*% t(H) %*% t(X) %*% W + beta * H %*% t(H) %*% Q %*% H + beta *( D %*% H- H %*% t(H)%*% D %*% H))  # 添加一个小的常数避免除零错误
    # 定义目标函数
    # tt[Ti+1]<-norm(X-W%*%t(H),"F")^2+alpha * sum(W * Marker)+beta*sum(diag(t(H) %*% L %*%H))+lambda*sum(H*A)
    # cha<-tt[Ti+1]-objective
    # objective<-tt[Ti+1]
    #KKT<-kkt[Ti+1]
    # H <- sweep(H, 1, rowSums(H) + 1e-8, "/")
    Ti<-Ti+1
    # print(cha)
  }
  print(Ti)
  print("Training Finished.")
  out<-list(W,H,objective,tt,Ti,cha)
  return(out)
}


genem2h <- function(dataset){
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
  mouse_genes <- rownames(dataset)
  
  # 映射
  human_mapped_genes <- mapping_df$human_symbol[
    match(mouse_genes, mapping_df$mouse_symbol)
  ]
  
  # 替换行名
  rownames(dataset) <- human_mapped_genes
  
  # 去掉没有匹配的人类基因
  dataset <- dataset[!is.na(rownames(dataset)), ]
  dataset <- dataset[rownames(dataset) != "", ]
  dataset <- dataset[!duplicated(rownames(dataset)), ]
  dataset
  }

metric_cal <- function(true_labels, annotation_labels){
  # 计算准确度
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
  
  # 打印每组的结果
  cat("Results for", ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")
  
  # 存储每组结果的指标
  list(
    true_labels = true_labels,
    pre_labels = annotation_labels,
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
}

# NNDSVD initialization for nonnegative matrix factorization
# X: nonnegative matrix (m x n)
# r: target rank
# variant: "nndsvd" (strict), "nndsvda" (fill zeros by the average of X),
#          "nndsvdar" (fill zeros with small random noise)
nndsvd_init <- function(X, r, variant = c("nndsvd", "nndsvda", "nndsvdar")) {
  variant <- match.arg(variant)
  stopifnot(all(X >= 0), r > 0, r <= min(dim(X)))
  
  m <- nrow(X); n <- ncol(X)
  
  # 用经济型 SVD；数据大时可替换成 irlba::irlba(X, r)
  s <- svd(X, nu = r, nv = r)
  U <- s$u[, 1:r, drop = FALSE]
  S <- s$d[1:r]
  V <- s$v[, 1:r, drop = FALSE]
  
  W <- matrix(0, m, r)
  H <- matrix(0, r, n)
  
  # 第一个分量：直接取正部
  W[, 1] <- sqrt(S[1]) * pmax(U[, 1], 0)
  H[1, ] <- sqrt(S[1]) * pmax(V[, 1], 0)
  
  # 为避免全零（极少见），如果某侧全 0，则用负部
  if (sum(W[,1]) == 0) W[,1] <- sqrt(S[1]) * pmax(-U[,1], 0)
  if (sum(H[1,]) == 0) H[1,] <- sqrt(S[1]) * pmax(-V[,1], 0)
  
  # 其余分量：按 Boutsidis & Gallopoulos (2008)
  for (k in 2:r) {
    u <- U[, k]; v <- V[, k]
    u_p <- pmax(u, 0); u_n <- pmax(-u, 0)
    v_p <- pmax(v, 0); v_n <- pmax(-v, 0)
    
    u_p_norm <- sqrt(sum(u_p^2)); u_n_norm <- sqrt(sum(u_n^2))
    v_p_norm <- sqrt(sum(v_p^2)); v_n_norm <- sqrt(sum(v_n^2))
    
    m_p <- u_p_norm * v_p_norm
    m_n <- u_n_norm * v_n_norm
    
    if (m_p >= m_n) {
      if (u_p_norm > 0) W[, k] <- u_p / u_p_norm else W[, k] <- 0
      if (v_p_norm > 0) H[k, ] <- v_p / v_p_norm else H[k, ] <- 0
      W[, k] <- W[, k] * sqrt(S[k] * m_p)
      H[k, ] <- H[k, ] * sqrt(S[k] * m_p)
    } else {
      if (u_n_norm > 0) W[, k] <- u_n / u_n_norm else W[, k] <- 0
      if (v_n_norm > 0) H[k, ] <- v_n / v_n_norm else H[k, ] <- 0
      W[, k] <- W[, k] * sqrt(S[k] * m_n)
      H[k, ] <- H[k, ] * sqrt(S[k] * m_n)
    }
  }
  
  # 处理零块：三种变体
  if (variant != "nndsvd") {
    zero_W <- (W == 0)
    zero_H = (H == 0)
    if (variant == "nndsvda") {
      avgX <- mean(X[X > 0])
      if (is.na(avgX) || avgX == 0) avgX <- 1e-6
      W[zero_W] <- avgX
      H[zero_H] <- avgX
    } else if (variant == "nndsvdar") {
      eps <- 1e-4 * mean(X[X > 0])
      if (is.na(eps) || eps == 0) eps <- 1e-4
      W[zero_W] <- runif(sum(zero_W), min = 0, max = eps)
      H[zero_H] <- runif(sum(zero_H), min = 0, max = eps)
    }
  }
  
  list(W = W, H = t(H))  # 注意：这里返回的 H 形状与你代码一致（n x r）
}

#基于单细胞转录组测序数据的细胞注释算法研究主代码
#X-基因表达矩阵，A标签矩阵，Marker标记基因矩阵，alpha-标记权重，beta-几何权重，lambda-标签权重，iteration-迭代次数，sigma-阈值
scANMF<-function(X,A,rank,Marker,alpha,beta,lambda,L,D,Q,iteration,sigma){
  # #initial W & H
  # kkk <- initial_NMF(X,rank,"nndsvd")
  # W <- kkk$W
  # H <- t(kkk$H)
  # print("Initialization Finished.")
  
  num_rows_W <- nrow(X)
  num_cols_W <- rank
  num_rows_H <- ncol(X)
  num_cols_H <- rank
  # 生成随机初始化的矩阵 U 和 Z
  W <- matrix(runif(num_rows_W * num_cols_W), nrow = num_rows_W, ncol = num_cols_W)
  H <- matrix(runif(num_rows_H * num_cols_H), nrow = num_rows_H, ncol = num_cols_H)
  # 
  # 打印生成的矩阵
  
  #计算目标函数
  objective<-norm(X-W%*%t(H),"F")^2+alpha * sum(W * Marker)+beta*sum(diag(t(H) %*% L %*%H))+lambda*sum(H*A)
  cha<-1
  Ti <- 0
  tt<-numeric(length = iteration)
  while (abs(cha)>sigma && Ti<=iteration){
    # print(objective)
    # print(norm(X-W%*%t(H),"F")^2)
    # print(alpha * sum(W * Marker))
    # print(beta*sum(diag(t(H) %*% L %*%H)))
    # print(lambda*sum(H*A))
    # 计算 X V
    XH <- X %*% H
    
    # 计算 UV^T V
    WHtH <- W %*% t(H) %*% H
    # 计算更新后的 U
    W= W * (2 * XH) / ((2 * WHtH + alpha*Marker))  # 添加一个小的常数避免除零错误
    # 计算更新后的 V
    H=H * (2*(t(X) %*% W  + 2*beta *  Q %*% H)) / (2*(H %*%t(W) %*% W +2* beta *  D %*% H)+lambda*A)
    #H = H * (t(X) %*% W + beta * Q %*% H) / ( H %*% t(H) %*% t(X) %*% W + beta * H %*% t(H) %*% Q %*% H + beta *( D %*% H- H %*% t(H)%*% D %*% H))  # 添加一个小的常数避免除零错误
    # 定义目标函数
    tt[Ti+1]<-norm(X-W%*%t(H),"F")^2+alpha * sum(W * Marker)+beta*sum(diag(t(H) %*% L %*%H))+lambda*sum(H*A)
    cha<-tt[Ti+1]-objective
    objective<-tt[Ti+1]
    #KKT<-kkt[Ti+1]
    Ti<-Ti+1
    print(cha)
  }
  print(Ti)
  print("Training Finished.")
  out<-list(W,H,objective,tt,Ti,cha)
  return(out)
}


mat_scaled <- ScaleData(seu, features = VariableFeatures(seu))@scale.data
mat_nonneg <- mat_scaled - matrixStats::colMins(mat_scaled)  # 按列平移

##最新梗概版本在这里
scANMF <- function(X, A, rank, Marker, alpha, beta, lambda, L, D, Q,
                   iteration, sigma, method_nmf = 'lee') {
  # 期望维度：
  # X: n×m,  W: n×r,  H: r×m
  # Marker: n×r,  A: r×m,  L/Q/D: m×m
  
  eps <- 1e-10

  
  
  # num_rows_W <- nrow(X)
  # num_cols_W <- rank
  # num_rows_H <- ncol(X)
  # num_cols_H <- rank
  # # 生成随机初始化的矩阵 U 和 Z
  # W <- matrix(runif(num_rows_W * rank), nrow = num_rows_W, ncol = rank)
  # H <- matrix(runif(rank * num_rows_H), nrow = rank, ncol = num_rows_H)
  
  library(NMF)
  res_nmf <- nmf(X, rank = rank, method_nmf, 'nndsvd', .opt = "vp")
  W <- basis(res_nmf)   # genes x k
  H <- coef(res_nmf)    # k x cells
  
  print(dim(A))
  print(dim(H))

  
  # 维度检查（出错就直接停）
  if (!(nrow(X) == nrow(W) && ncol(X) == ncol(H) && ncol(W) == nrow(H))) {
    stop("Dimension mismatch: ensure X(n×m), W(n×r), H(r×m).")
  }
  if (!all(dim(Marker) == dim(W))) stop("Marker must be n×r (same as W).")
  if (!all(dim(A) == dim(H)))      stop("A must be r×m (same as H).")
  if (!(nrow(L) == ncol(L) && nrow(L) == ncol(H))) stop("L must be m×m.")
  if (!(nrow(Q) == ncol(Q) && nrow(Q) == ncol(H))) stop("Q must be m×m.")
  if (!(nrow(D) == ncol(D) && nrow(D) == ncol(H))) stop("D must be m×m.")
  
  message("Initialization Finished.")
  
  ## 2) 目标函数（安全计算）
  frob2 <- function(M) sum(M * M)
  
  objective <- frob2(X - W %*% H) +
    alpha  * sum(W * Marker) +
    beta   * sum(diag(H %*% L %*% t(H))) +
    lambda * sum(H * A)
  
  cha <- 1
  Ti  <- 0
  tt  <- numeric(length = iteration)
  
  ## 3) 迭代
  while (abs(cha) > sigma && Ti < iteration) {
    # --- 3.1 更新 W ---
    # Numerator: X %*% t(H)  -> n×r
    # Denominator: W %*% H %*% t(H) + (alpha*Marker)/2 -> n×r
    XHt  <- X %*% t(H)
    WHHt <- W %*% H %*% t(H)
    
    W_num <- XHt
    W_den <- WHHt + (alpha * Marker) / 2 + eps
    
    W <- W * (W_num / W_den)
    # 数值稳定：去除非有限，截断为非负
    W[!is.finite(W)] <- 0
    W <- pmax(W, eps)
    
    # --- 3.2 更新 H ---
    # Numerator: t(W) %*% X + 2*beta * (H %*% Q)  -> r×m
    # Denominator: t(W) %*% W %*% H + 2*beta * (H %*% D) + lambda*A -> r×m
    tW    <- t(W)
    tWX   <- tW %*% X
    tWWH  <- tW %*% W %*% H
    
    H_num <- tWX + beta * (H %*% Q)
    H_den <- tWWH + beta * (H %*% D) + lambda * A/2 + eps
    
    H <- H * (H_num / H_den)
    H[!is.finite(H)] <- 0
    H <- pmax(H, eps)
    
    # --- 3.3 重新计算目标函数与收敛量 ---
    term_recon  <- frob2(X - W %*% H)
    term_marker <- alpha  * sum(W * Marker)
    term_graph  <- beta   * sum(diag(H %*% L %*% t(H)))
    term_mask   <- lambda * sum(H * A)
    
    tt[Ti + 1] <- term_recon + term_marker + term_graph + term_mask
    
    cha <- tt[Ti + 1] - objective
    objective <- tt[Ti + 1]
    Ti <- Ti + 1
    
    # 若出现非有限值，提前退出以避免 NaN 传播
    if (!is.finite(objective) || is.na(objective)) {
      warning(sprintf("Non-finite objective at iter %d; stopping.", Ti))
      break
    }
    # message(sprintf("iter=%d, Δ=%.4e, obj=%.6e", Ti, cha, objective))
  }
  
  message("Training Finished.")
  
  # 返回具名列表
  return(list(
    W = W, H = t(H),
    objective = objective,
    objective_trace = tt[seq_len(Ti)],
    iters = Ti,
    last_delta = cha
  ))
}
