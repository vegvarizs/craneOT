# R/kernel.R
#' Build sparse kNN kernel (exponential)
#' Build sparse k-NN Gibbs kernel
#' @param coords numeric matrix or data frame (eg. \code{lon}, \code{lat} columns)
#' @param k      integer; number of neighbours
#' @param eps    numeric; distance scale
#'
#' @return lista: \code{K} (dgCMatrix, \eqn{n \times n}) és \code{D} (distance matrix)
#' @export
build_sparse_kernel <- function(coords, k = 60, eps = 50) {
  stopifnot(ncol(coords) >= 2)
  XY <- as.matrix(coords[,1:2, drop=FALSE])
  n  <- nrow(XY); if (n < 2) stop("Need >=2 sites")
  D  <- as.matrix(dist(XY))
  diag(D) <- 0
  nbrs <- lapply(seq_len(n), function(i){
    ord <- order(D[i,])
    ord[ord != i][1:min(k, n-1)]
  })
  edges <- do.call(rbind, lapply(seq_len(n), function(i){
    j <- nbrs[[i]]
    if (length(j)) rbind(cbind(i,j), cbind(j,i))
  }))
  edges <- unique(edges)
  w <- exp(- D[edges] / eps)
  K <- Matrix::sparseMatrix(
    i = c(edges[,1], seq_len(n)),
    j = c(edges[,2], seq_len(n)),
    x = c(w, rep(1, n)),
    dims = c(n, n)
  )
  list(K = methods::as(K, "dgCMatrix"), D = D)
}
