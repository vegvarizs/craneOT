# R/sinkhorn.R
#' Unbalanced entropic Sinkhorn on sparse kernel
#' @param a,b numeric nonnegative marginals (length n)
#' @param K dgCMatrix kernel
#' @param eps,rho numeric, regularisation and marginal penalty
#' @param maxiter,tol control
#' @return dgCMatrix transport plan
#' @export
sinkhorn_uot <- function(a, b, K, eps=60, rho=500, maxiter=2000, tol=1e-8){
  a <- pmax(as.numeric(a),0); b <- pmax(as.numeric(b),0)
  n <- length(a); stopifnot(n==nrow(K), n==ncol(K))
  tau <- rho/(rho+eps); tiny <- .Machine$double.eps
  u <- rep(1,n); v <- rep(1,n); K <- methods::as(K,"dgCMatrix")
  for (it in seq_len(maxiter)) {
    Kv  <- as.numeric(K %*% v);  u_new <- (a/pmax(Kv,tiny))^tau
    Ktu <- as.numeric(Matrix::t(K) %*% u_new); v_new <- (b/pmax(Ktu,tiny))^tau
    if (max(abs(u_new-u), abs(v_new-v)) < tol) { u<-u_new; v<-v_new; break }
    u<-u_new; v<-v_new
  }
  F <- K; F@x <- F@x * rep(v, diff(F@p)); F <- Matrix::Diagonal(x=u) %*% F
  methods::as(F, "dgCMatrix")
}
