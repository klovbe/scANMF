scANMF<-function(X,A,rank,Marker,alpha,beta,lambda,L,D,Q,iteration,sigma,init='runif'){
  ##initial W & H
  set.seed(123)
  num_rows_W <- nrow(X)
  num_cols_W <- rank
  num_rows_H <- ncol(X)
  num_cols_H <- rank
  if(init == 'runif'){
    W <- matrix(runif(num_rows_W * num_cols_W), nrow = num_rows_W, ncol = num_cols_W)
    H <- matrix(runif(num_rows_H * num_cols_H), nrow = num_rows_H, ncol = num_cols_H)
  } 

  else{
    init <- nndsvd_init(X, rank, variant = "nndsvda")  # or "nndsvd"/"nndsvdar"
    W <- init$W
    H <- init$H  
  }
  
  
  MM <- alpha*Marker
  zero_rows <- rowSums(MM == 0) == ncol(MM)
  MM[zero_rows, ] <- 50
  AA <- lambda*A

  
  
  # objective<-norm(X-W%*%t(H),"F")^2+alpha * sum(W * Marker)+beta*sum(diag(t(H) %*% L %*%H))+lambda*sum(H*A)
  cha<-1
  Ti <-0
  objective<-0
  tt<-numeric(length = iteration)
  while (abs(cha)>sigma && Ti<=iteration){
    XH <- X %*% H
    WHtH <- W %*% t(H) %*% H
    W= W * (2 * XH) / ((2 * WHtH + MM)) 
    H=H * (2*(t(X) %*% W  + 2*beta *  Q %*% H)) / (2*(H %*%t(W) %*% W +2* beta *  D %*% H)+AA)
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
  
  orthologs <- homologene::homologeneData
  colnames(orthologs)
  orthologs_hm <- subset(orthologs, Taxonomy %in% c(9606, 10090))

  orthologs_hm <- reshape(
    orthologs_hm[, c("HID", "Gene.Symbol", "Taxonomy")],
    timevar = "Taxonomy",
    idvar = "HID",
    direction = "wide"
  )
  

  colnames(orthologs_hm) <- c("HomoloGeneID", "human_symbol", "mouse_symbol")
  mapping_df <- orthologs_hm[, c("mouse_symbol", "human_symbol")]

  mouse_genes <- rownames(dataset)
  human_mapped_genes <- mapping_df$human_symbol[
    match(mouse_genes, mapping_df$mouse_symbol)
  ]

  rownames(dataset) <- human_mapped_genes

  dataset <- dataset[!is.na(rownames(dataset)), ]
  dataset <- dataset[rownames(dataset) != "", ]
  dataset <- dataset[!duplicated(rownames(dataset)), ]
  dataset
  }

metric_cal <- function(true_labels, annotation_labels){
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
  
  cat("Results for", ":\n")
  cat("Accuracy:", accuracy, "\n")
  cat("Macro F1:", mean_f1_score, "\n")
  cat("Weighted F1:", weighted_f1, "\n\n")

  list(
    true_labels = true_labels,
    pre_labels = annotation_labels,
    accuracy = accuracy,
    macro_f1 = mean_f1_score,
    weighted_f1 = weighted_f1,
    f1_score_per_class = f1_score
  )
}
