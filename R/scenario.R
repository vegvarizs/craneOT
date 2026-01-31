#' Rerun selected weeks under a modified kernel
#'
#' Recomputes unbalanced-OT flows for chosen \eqn{t \to t+1} transitions
#' using a modified sparse kernel.
#'
#' @param A Numeric matrix of dimension \eqn{T \times n}; rows are weeks,
#'   columns are roosts (observed counts).
#' @param weeks_idx Integer vector of week indices \eqn{t} for which both
#'   \eqn{t} and \eqn{t+1} exist.
#' @param K_mod \code{dgCMatrix} sparse kernel of size \eqn{n \times n}.
#' @param eps Positive numeric; entropic regularisation parameter \eqn{\varepsilon}.
#' @param rho Positive numeric; marginal relaxation penalty \eqn{\rho}.
#' @param sinkhorn_fun Function with signature
#'   \code{function(a_t, b_t, K, eps, rho, ...)}; defaults to \code{\link{sinkhorn_uot}}.
#' @param maxiter Integer; maximum number of iterations.
#' @param tol Numeric; convergence tolerance.
#'
#' @return A list of \eqn{n \times n} (usually sparse) flow matrices \eqn{F_t},
#'   one for each requested transition.
#'
#' @details For each \eqn{t \in \mathrm{weeks\_idx}}, the function computes
#'   \eqn{a_t = A[t,\,]} and \eqn{b_t = A[t+1,\,]}, then calls
#'   \code{sinkhorn_fun(a_t, b_t, K_mod, eps, rho, maxiter = maxiter, tol = tol)}.
#'
#' @examples
#' # out <- rerun_weeks(A, weeks_idx = c(10, 11), K_mod, eps = 60, rho = 30)
#' @export
rerun_weeks <- function(A, weeks_idx, K_mod, eps, rho,
                        sinkhorn_fun = sinkhorn_uot,
                        maxiter = 1000, tol = 1e-6) {
  n <- ncol(A)
  stopifnot(nrow(K_mod) == n, ncol(K_mod) == n)
  out <- vector("list", length(weeks_idx))
  for (k in seq_along(weeks_idx)) {
    t <- weeks_idx[k]
    a_t <- A[t, ]
    b_t <- A[t + 1, ]
    out[[k]] <- sinkhorn_fun(a_t, b_t, K_mod, eps = eps, rho = rho,
                             maxiter = maxiter, tol = tol)
  }
  out
}

