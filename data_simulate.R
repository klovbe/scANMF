model_sample <- function(K, Kn_vec, Ndiff, Nsame, logMean, logSd, ZeroRate, sigmahetero,sigmahomo,type, a, b, error_rate) {
  K <- K
  N <- sum(Kn_vec)  # Total number of cells
  P <- Ndiff + Nsame  # Total number of genes
  
  # Generate the labels with an uneven distribution based on Kn_vec
  Label <- unlist(mapply(function(k, n) rep(k, n), 1:K, Kn_vec))  # Uneven labels
  
  Z <- matrix(0, N, K)  # Initialize Z matrix
  start <- 1
  for (i in 1:K) {
    end <- start + Kn_vec[i] - 1
    Z[start:end, i] <- 1  # Assign Z values for each cell type
    start <- end + 1
  }
  
  # Generate expression values for genes
  Esame <- rnorm(Nsame, logMean, logSd)
  Esame <- Esame * rbinom(Nsame, 1, 1 - ZeroRate)  # Simulate true zeros
  Esame <- rep(Esame, K)
  dim(Esame) <- c(Nsame, K)
  
  if (type == "cluster") {
    Ediff <- matrix(rnorm(K * Ndiff, logMean, logSd), nrow = K, ncol = Ndiff)
  } else if (type == "DE") {
    Ediff <- rnorm(Ndiff, logMean, logSd)
    Ediff <- matrix(Ediff, nrow = K, ncol = Ndiff, byrow = TRUE)
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] <- Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] * 2
    }
  } else if (type == "marker") {
    Ediff <- matrix(rnorm(K * Ndiff, logMean - 0.5, logSd), nrow = K, ncol = Ndiff)
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      Ediff[k, ((k - 1) * Ndiff.u + 1):(k * Ndiff.u)] <- rnorm(Ndiff.u, logMean, logSd) * 1.5
    }
  }
  Ediff <- Ediff * matrix(rbinom(K * Ndiff, 1, 1 - ZeroRate), K, Ndiff)
  
  Te <- Z %*% cbind(Ediff, t(Esame))
  Te[Te > 0] <- Te[Te > 0] + rnorm(sum(Te > 0), 0, Te[Te > 0] * sigmahetero + sigmahomo)
  Te[Te <= 0] <- rnorm(sum(Te <= 0), 0, 1) - 1.5
  Te[Te <= 0] <- 0
  
  De <- t(Te)
  
  # Marker gene matrix
  markerGenes <- matrix(1, nrow = K, ncol = P)  # Initialize as 1
  if (type == "marker") {
    Ndiff.u <- round(Ndiff / K)
    for (k in 1:K) {
      for (j in 1:Ndiff.u) {
        markerGenes[k, ((k - 1) * Ndiff.u + j)] <- 0  # Set marker gene positions to 0
      }
    }
  }
  
  # Generate sample, type, and gene names
  sample_names <- paste("Sample", 1:N, sep = "")
  type_names <- paste("Type", 1:K, sep = "")
  gene_names <- paste("Gene", 1:P, sep = "")
  
  rownames(markerGenes) <- type_names
  colnames(markerGenes) <- gene_names
  markerGenes <- t(markerGenes)
  rows_all_ones <- apply(markerGenes, 1, function(x) all(x == 1))
  markerGenes[rows_all_ones, ] <- 0
  
  rownames(Te) <- sample_names
  colnames(Te) <- gene_names
  
  # Marker gene information
  marker_gene_info <- list()
  for (i in 1:P) {
    marker_in_types <- which(markerGenes[i, ] == 0)
    if (length(marker_in_types) > 0) {
      marker_gene_info[[gene_names[i]]] <- type_names[marker_in_types]
    }
  }
  
  # Transpose gene expression matrix
  Te <- t(Te)
  
  # Generate correct label matrix
  num_cell_types <- length(unique(Label))
  label_matrix_correct <- matrix(0, nrow = N, ncol = num_cell_types)
  for (i in 1:num_cell_types) {
    label_matrix_correct[Label == i, i] <- 1
  }
  
  # Generate label matrix with errors
  label_matrix_with_errors <- label_matrix_correct
  percentage <- a
  for (i in 1:num_cell_types) {
    cell_indices <- which(Label == i)
    num_to_keep <- ceiling(percentage * length(cell_indices))
    keep_indices <- sample(cell_indices, num_to_keep)
    drop_indices <- setdiff(cell_indices, keep_indices)
    label_matrix_with_errors[drop_indices, ] <- 0
    label_matrix_with_errors[keep_indices, i] <- 1
  }
  
  for (i in 1:num_cell_types) {
    cell_indices <- which(label_matrix_with_errors[, i] == 1)
    num_errors <- ceiling(b * length(cell_indices))
    error_indices <- sample(cell_indices, num_errors)
    for (idx in error_indices) {
      other_types <- setdiff(1:num_cell_types, i)
      label_matrix_with_errors[idx, i] <- 0
      label_matrix_with_errors[idx, sample(other_types, 1)] <- 1
    }
  }
  
  celltype_matrix <- label_matrix_with_errors
  for (i in 1:N) {
    if (all(celltype_matrix[i, ] == 0)) {
      celltype_matrix[i, ] <- label_matrix_correct[i, ]
    }
  }
  
  # Generate marker genes with errors
  generate_marker_with_errors <- function(n_genes, n_cell_types, error_rate) {
    total_elements <- n_genes * n_cell_types
    num_errors <- round(total_elements * error_rate)
    
    error_indices <- sample(1:total_elements, num_errors)
    markerGenes_error <- markerGenes
    for (index in error_indices) {
      row_index <- ceiling(index / n_cell_types)
      col_index <- index %% n_cell_types
      if (col_index == 0) col_index <- n_cell_types
      markerGenes_error[row_index, col_index] <- 1 - markerGenes_error[row_index, col_index]
    }
    return(markerGenes_error)
  }
  
  markerGenes_error <- generate_marker_with_errors(P, K, error_rate)
  
  # Return simulation data
  sData <- list(K = K, N = N, P = P, label = Label, te = Te, de = De, markerGenes_error = markerGenes_error,
                labelMatrixCorrect = label_matrix_correct, labelMatrixWithErrors = label_matrix_with_errors,
                celltype_matrix = celltype_matrix, marker_gene_info = marker_gene_info)
  return(sData)
}

  
