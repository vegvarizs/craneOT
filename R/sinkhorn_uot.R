
#' Sparse unbalanced Sinkhorn OT
#' @param a,b numeric vectors (length n)
#' @param K   sparse kernel (dgCMatrix, n x n)
#' @param eps,rho positive scalars
#' @return sparse transport plan (dgCMatrix)
#' @export
sinkhorn_uot <- function(a, b, K, eps = 60, rho = 30, maxiter = 1500, tol = 1e-6){
  stopifnot(inherits(K, "dgCMatrix"))
  a <- pmax(as.numeric(a), 0); b <- pmax(as.numeric(b), 0)
  n <- length(a); if (n != nrow(K) || n != ncol(K)) stop("K dimension mismatch")
  tau <- rho / (rho + eps); tiny <- .Machine$double.eps
  u <- v <- rep(1, n)
  for (it in seq_len(maxiter)) {
    Kv  <- as.numeric(K %*% v);  u_new <- (a / pmax(Kv, tiny))^tau
    Ktu <- as.numeric(Matrix::t(K) %*% u_new); v_new <- (b / pmax(Ktu, tiny))^tau
    if (max(abs(u_new-u), abs(v_new-v)) < tol) { u <- u_new; v <- v_new; break }
    u <- u_new; v <- v_new
  }
  Matrix::Diagonal(x = u) %*% K %*% Matrix::Diagonal(x = v)
}

