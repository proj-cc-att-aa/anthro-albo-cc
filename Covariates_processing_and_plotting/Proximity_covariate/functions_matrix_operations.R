# ================== Matrix Operations ==================

packages <- c("Matrix")
lapply(packages, library, character.only = TRUE)




# eigen_decomp.
eigen_decomp <- function(A) {
  
  n <- nrow(A)
  rank <- rankMatrix(A)
  E <- eigen(A)
  determinant_log <- sum(log(E$values[1:rank]))
  
  colnames(E$vectors) <- sapply(1:n, function(x) paste0("X", x))
  
  return(list(n = n,
              rank = rank,
              determinant_log = determinant_log,
              E = E))
}



#plot_eigenvector.
plot_eigenvector <- function(Q_eigen, map, i, title = "") {
  eigen_vec <- Q_eigen$E$vectors[,i]
  eigen_val <- Q_eigen$E$values[i]
  descr <- paste0("eig.vec.: ", i,
                  ", eig.val.: ", eigen_val %>%
                    sprintf("%.2f", .))
  
  map_eigen <- cbind(map, eigen_vec = eigen_vec)
  
  field = "eigen_vec"
  if (title == "") {
    title = descr
  } else {
    title = paste0(title, " - ", descr)
  }
  
  if ("geometry" %in% names(map)) {
    p <- plot_sf(map_eigen, field, title = title) %>% print()
  } else if (all(is.element(c("longitude", "latitude"), names(map)))) {
    p <- plot_scatter_lines(map_eigen, "longitude", "latitude", field,
                            plot_lines=FALSE, title=title)
  } else {
    stop("eigen-plot error")
  }
  
  return(p)
}

#  matrix_regular_to_eig_basis.
matrix_regular_to_eig_basis <- function(M, Q_eigen) {
  M_E <- t(Q_eigen$E$vectors) %*% M %*% Q_eigen$E$vectors
  
  return(M_E)
}

# matrix_eig_to_regular_basis.
matrix_eig_to_regular_basis <- function(M_E, Q_eigen) {
  M <- Q_eigen$E$vectors %*% M_E %*% t(Q_eigen$E$vectors)
  
  return(M)
}

# plot_matrix.
plot_matrix <- function(x, reorder = FALSE,
                        reorder_rows = FALSE, reorder_cols = FALSE,
                        title = "", palette_idx = 2, col_breaks = NULL,
                        title_font_sc = 0.75,
                        label_font_size = NULL) {
  
  x <- as.matrix(x)
  
  cols <- color_palette(palette_idx,)[["col_palette"]]
  pal <- colorRampPalette(cols)(50)
  
  title_font_par <- par()$cex.main
  title_font_matrix <- title_font_par * title_font_sc
  par(cex.main = title_font_matrix)
  
  if (is.null(label_font_size)) {
    label_font_size_rows <- 0.2 + 1/log10(nrow(x))
    label_font_size_cols <- 0.2 + 1/log10(ncol(x))
  } else {
    label_font_size_rows <- label_font_size
    label_font_size_cols <- label_font_size
  }
  
  if (reorder) {
    if (is.null(col_breaks)) {
      x_hm <- heatmap.2(x, scale="none",
                        main=title,
                        xlab="Columns", ylab="Rows",
                        col=pal, tracecol="#303030", trace="none",
                        keysize = 1.5, margins=c(4, 4), symbreaks=FALSE,
                        cexRow = label_font_size_rows, cexCol = label_font_size_cols)
    } else {
      x_hm <- heatmap.2(x, scale="none",
                        Rowv=FALSE, Colv=FALSE, dendrogram="none",
                        main=title,
                        xlab="Columns", ylab="Rows",
                        col=pal, tracecol="#303030", trace="none",
                        breaks = col_breaks,
                        keysize = 1.5, margins=c(4, 4),
                        cexRow = label_font_size_rows, cexCol = label_font_size_cols)
    }
  } else {
    
    if (reorder_rows) {
      Rowv = TRUE
      Colv = FALSE
    } else if (reorder_cols) {
      Rowv = FALSE
      Colv = TRUE
    } else {
      Rowv = FALSE
      Colv = FALSE
    }
    
    if (is.null(col_breaks)) {
      x_hm <- heatmap.2(x, scale="none",
                        Rowv=Rowv, Colv=Colv, dendrogram="none",
                        main=title,
                        xlab="Columns", ylab="Rows",
                        col=pal, tracecol="#303030", trace="none",
                        keysize = 1.5, margins=c(4, 4),
                        cexRow = label_font_size_rows, cexCol = label_font_size_cols)
    } else {
      x_hm <- heatmap.2(x, scale="none",
                        Rowv=Rowv, Colv=Colv, dendrogram="none",
                        main=title,
                        xlab="Columns", ylab="Rows",
                        col=pal, tracecol="#303030", trace="none",
                        breaks = col_breaks,
                        keysize = 1.5, margins=c(4, 4),
                        cexRow = label_font_size_rows, cexCol = label_font_size_cols)
    }
  }
  
  par(cex.main = title_font_par)
  
  return(x_hm)
}



# simulate_from_eig_decomp.
simulate_from_eig_decomp <- function(eig_vectors, eig_values) {
  
  n <- length(eig_values)
  
  z_weighted <- rnorm(n,mean=0,sd=1) * eig_values
  x <- eig_vectors %*% z_weighted
  
  return(list(x=x, z_weighted=z_weighted))
}




#  simulate_from_chol_decomp.
simulate_from_chol_decomp <- function(Q_chol=NULL, C_chol=NULL) {
  
  if (!is.null(Q_chol)) {
    
    n <- nrow(Q_chol)
    
    x <- as.vector(solve(Q_chol, rnorm(n,mean=0,sd=1)))
  } else if (!is.null(C_chol)) {
    
    n <- nrow(C_chol)
    
    x <- as.vector(t(C_chol) %*% rnorm(n,mean=0,sd=1))
  }
  
  return(x)
}


#  simulate_with_inla.
simulate_with_inla <- function(Q) {
  
  x <- as.vector(inla.qsample(n=1, Q))
  
  return(x)
}



# generalized_inverse_through_eig_decomp.
generalized_inverse_through_eig_decomp <- function(eig_vectors, eig_values) {
  
  n <- length(eig_values)
  
  d <- rep(0, n)
  
  eig_values_standardized <- eig_values / sd(eig_values)
  if (any(abs(eig_values_standardized) < 1e-10)) {
    
    i_zero <- which(abs(eig_values_standardized) < 1e-10)
    
    d[i_zero] <- 0
    d[-i_zero] <- 1 / eig_values[-i_zero]
  } else {
    
    d <- 1 / eig_values
  }
  
  generalized_inverse <- eig_vectors %*% diag(d) %*% t(eig_vectors)
  
  return(generalized_inverse)
}
