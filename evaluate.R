library(Matrix)
library(gplots)

# reconstruction error of each pattern
pat_recon_err <- function(X, Y, ground_truth) {
  "
  This function computes the mean reconstruction error of patterns
  Input:
    X, Y:           factor matrices got by a BMF method
    ground_truth:   ground truth pattern matrix
  Output:
    error:          pattern reconstruction error
  "
  mat_dim <- dim(X)[[1]]
  patterns <- as.matrix(X %&% Y * 1)
  error_nums <- sum(abs(patterns - ground_truth))
  error <- error_nums / (mat_dim*mat_dim)
  return(list(error_nums, error))
}

# estimation error of pr and pc
bg_accuracy <- function(p_rows, p_cols, true_prows, true_pcols) {
  "
  This function computes errors of the estimated background probs
  Input:
    p_rows, p_cols:         probs estimated by the denoising algorithm
    true_prows, true_pcols: ground truth probs
  "
  error_rows <- mean(abs(p_rows - true_prows))
  error_cols <- mean(abs(p_cols - true_pcols))
  return(list(error_rows, error_cols))
}


denoise_eval <- function(X_pattern, X_bg, mat_after_denoise) {
  # metrics about signal and noise
  bg_pat_overlap <- X_pattern&X_bg
  X_bg_bool <- X_bg > 0
  X_pat_bool <- X_pattern > 0
  pure_bg <- (0+X_bg_bool)- (0+bg_pat_overlap)
  pure_pat <- (0+X_pat_bool) - (0+bg_pat_overlap)
  bg_strength <- sum(pure_bg)
  pat_strength <- sum(X_pattern)
  sn_ratio <- pat_strength / bg_strength  # original signal noise ratio
  
  pat_after_denoise <- X_pattern&mat_after_denoise
  bg_after_denoise <- pure_bg&mat_after_denoise
  overlap_after_denoise <- bg_pat_overlap&mat_after_denoise
  pure_pat_left <- pure_pat&mat_after_denoise
  bg_strength2 <- sum(bg_after_denoise)
  pat_strength2 <- sum(pat_after_denoise)
  denoised_sn_ratio <- pat_strength2 / bg_strength2  # signal noise ratio of denoised data
                                                     
  retained_pat <- pat_strength2 / pat_strength  # percentage of patterns retained
  
  # metrics about denoise accuracy
  TP <- (bg_strength - bg_strength2) / bg_strength
  # FP <- (pat_strength - pat_strength2) / pat_strength
  FP1 <- (sum(bg_pat_overlap)-sum(overlap_after_denoise)) / sum(bg_pat_overlap)
  FP2 <- (sum(pure_pat)-sum(pure_pat_left)) / sum(pure_pat)
  FN <- bg_strength2 / bg_strength
  TN <- pat_strength2 / pat_strength
  precision <- TP / (TP+FP1+FP2)
  precision2 <- (TP + FP1) / (TP+FP1+FP2)
  recall <- TP / (TP+FN)
  
  res = list(sn_ratio, denoised_sn_ratio, retained_pat, precision, precision2, recall)
  names(res) = c("sn_ratio", "denoised_sn_ratio", "retained_pat", "precision", "precision2", "recall")
  return(res)
}

jaccard_index <- function(X_data, X_pattern) {
  common_matrix <- X_data&X_pattern
  jaccard <- sum(common_matrix) / 
    (sum(X_data) + sum(X_pattern) - sum(common_matrix))
  return(jaccard)
}


myplot <- function(data, patterns) {
  colors = c(-100:100) / 100
  # my_palette <- colorRampPalette(c("#E41A1C", "white", "#377EB8"))(n = 200)
  my_palette <- colorRampPalette(c("brown1", "white", "dodgerblue2"))(n = 200)
  xx0 = patterns
  xx2 = data
  xx2[which((xx2 == 1) & (xx0 == 1))] <- -1
  # h <- heatmap.2(xx2, Rowv = T, Colv = T, scale = "none", main = "", 
  #                col = my_palette, breaks = colors, density.info = "none",
  #                dendrogram = "none", trace = "none", margin = c(10, 10),
  #                cexRow = 0.5,cexCol = 0.5, key = FALSE)
  h <- heatmap.2(xx2, col = my_palette, breaks = colors, density.info = "none",
                 dendrogram = "none", trace = "none", margin = c(10, 10),
                 cexRow = 0.5,cexCol = 0.5, key = FALSE)
}


myplot2 <- function(data) {
  my_palette <- colorRampPalette(c("white", "blue"))(100)
  h <- heatmap.2(data, trace = "none", col = my_palette, Rowv = FALSE, 
                 Colv = FALSE, dendrogram = "none", key = FALSE, 
                 margins = c(10,10))
}


