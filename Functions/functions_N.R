# ==================Computation of Q ==================

library("gplots")


Q_improper_by_v <- function(C, v, n_iter=100) {  
  
  alpha <- as.numeric(t(v) %*% v)
  
  C_diff_v <- C %*% (diag(n_regions) - v %*% (t(v) / alpha))
  
  print(t(v) %*% C_diff_v %*% v)
  
  C_diff_v_eigen <- eigen_decomp(C_diff_v)
  Q_diff_v <- 
    generalized_inverse_through_eig_decomp(C_diff_v_eigen$E$vectors, 
                                           C_diff_v_eigen$E$values)
  
  # # Q_v_eigen <- eigen_decomp(Q_diff_v)
  # 
  # # v_eigen <- t(Q_v_eigen$E$vectors) %*% v
  # # 
  # # plot(v_eigen)
  # # plot((Q_v_eigen$E$values)^0.1)
  # 
  # if (n_iter > 0) {
  #   for (j in 1:n_iter) {
  #     Q_v_eigen <- eigen_decomp(Q_diff_v)
  #     
  #     Q_diff_v_from_eig <- Q_v_eigen$E$vectors %*% 
  #       diag(Q_v_eigen$E$values) %*% 
  #       t(Q_v_eigen$E$vectors)
  #     
  #     t(v) %*% Q_diff_v_from_eig %*% v
  #     
  #     Q_diff_v_round2 <- Q_diff_v_from_eig %*% (diag(n_regions) - v %*% (t(v) / alpha))
  #     
  #     t(v) %*% Q_diff_v_round2 %*% v
  #     
  #     Q_v_eigen <- eigen_decomp(Q_diff_v_round2)
  #     
  #     v_eigen <- t(Q_v_eigen$E$vectors) %*% v
  #     
  #     par(mfrow=c(1,2))
  #     plot(v_eigen)
  #     plot((Q_v_eigen$E$values)^0.1)
  #     
  #     Q_diff_v <- Q_diff_v_round2
  #   }
  # }
  plot_matrix(C_diff_v - C, reorder=FALSE)
  # plot_matrix(t(t(Q_diff_v / diag(Q_diff_v)) / diag(Q_diff_v)), reorder=FALSE)
  # plot_matrix(t(t(Q / diag(Q)) / diag(Q)), reorder=FALSE)
  
  # par(mfrow=c(1,1))
  
  return(Q_diff_v)
}