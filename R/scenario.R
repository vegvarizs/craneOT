# R/scenario.R
#' Rerun selected weeks under modified kernel
#' @export
rerun_weeks <- function(A, weeks_idx, K_mod, eps, rho){
  n_out <- length(weeks_idx)
  out <- vector("list", n_out)
  for (i in seq_along(weeks_idx)) {
    t <- weeks_idx[i]
    out[[i]] <- sinkhorn_uot(A[,t], A[,t+1], K_mod, eps=eps, rho=rho)
  }
  out
}
