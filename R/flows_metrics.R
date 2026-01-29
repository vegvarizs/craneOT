# R/flows_metrics.R
#' Aggregate flows and compute hub scores & mean distance
#' @param flows list of n x n (dgCMatrix) weekly flows
#' @param Dmat optional n x n distances for mean distance
#' @param A optional n x T counts to scale hub score
#' @export
flows_summary <- function(flows, Dmat=NULL, A=NULL){
  Psum <- Reduce(`+`, flows)
  inflow  <- as.numeric(Matrix::colSums(Psum))
  outflow <- as.numeric(Matrix::rowSums(Psum))
  hub     <- inflow + outflow
  scale_tot <- if (!is.null(A)) sum(A) else sum(hub)
  hub_score <- hub / max(scale_tot, 1)
  md <- NA_real_
  if (!is.null(Dmat)) {
    num <- sum(Psum * Dmat); den <- sum(Psum); md <- ifelse(den>0, num/den, NA_real_)
  }
  list(Psum=Psum, inflow=inflow, outflow=outflow, hub_score=hub_score, mean_distance=md)
}

#' External flux for one transition
#' @export
external_flux <- function(Ft, a_t, b_t){
  ext_out <- sum(pmax(a_t - as.numeric(Matrix::rowSums(Ft)), 0))
  ext_in  <- sum(pmax(b_t - as.numeric(Matrix::colSums(Ft)), 0))
  c(ext_out=ext_out, ext_in=ext_in)
}
