#' Flow-derived edge and hub metrics
#' Summaries from a list of weekly flow matrices
#'
#' Computes cumulative flows (\eqn{P_{\mathrm{sum}}=\sum_t F_t}),
#' in-/outflows, a hub score, and—optionally—the global mean transport
#' distance if a distance matrix is provided.
#'
#' @param flows List of weekly \eqn{n \times n} numeric or sparse
#'   (\code{dgCMatrix}) flow matrices \eqn{F_t}.
#' @param Dmat Optional \eqn{n \times n} distance matrix \eqn{d_{ij}}.
#'   If supplied, the mean transport distance is
#'   \eqn{\sum_{ij} P_{\mathrm{sum},ij} d_{ij} / \sum_{ij} P_{\mathrm{sum},ij}}.
#' @param A Optional \eqn{T \times n} count matrix. If supplied, the
#'   \code{hub_score} is scaled by \code{sum(A)}; otherwise by
#'   \code{sum(inflow + outflow)}.
#'
#' @return A list with:
#' \describe{
#'   \item{Psum}{\eqn{n \times n} sparse matrix; cumulative flow \eqn{\sum_t F_t}.}
#'   \item{inflow}{Numeric length-\eqn{n}; \eqn{\sum_j P_{\mathrm{sum},ji}}.}
#'   \item{outflow}{Numeric length-\eqn{n}; \eqn{\sum_j P_{\mathrm{sum},ij}}.}
#'   \item{hub_score}{Numeric length-\eqn{n};
#'         \eqn{(inflow+outflow)/\max(1,\mathrm{scale\_tot})}.}
#'   \item{mean_distance}{Scalar or \code{NA};
#'         \eqn{\sum_{ij} P_{\mathrm{sum},ij} d_{ij} / \sum_{ij} P_{\mathrm{sum},ij}}
#'         if \code{Dmat} is supplied.}
#' }
#'
#' @examples
#' # res <- flows_summary(flows)
#' # str(res$Psum)
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

#' External import/export for one transition
#'
#' Computes network-external exports and imports for a single
#' \eqn{t \to t+1} transition from a reconstructed flow matrix and
#' observed margins.
#'
#' @param Ft \eqn{n \times n} numeric or sparse (\code{dgCMatrix}) flow matrix
#'   for the \eqn{t \to t+1} transition.
#' @param a_t Numeric length-\eqn{n}; observed totals at week \eqn{t}.
#' @param b_t Numeric length-\eqn{n}; observed totals at week \eqn{t+1}.
#'
#' @return A list with two scalars:
#' \describe{
#'   \item{ext_out}{\eqn{\sum_i \max\{0,\, a_{i,t} - (F_t \mathbf 1)_i\}}.}
#'   \item{ext_in }{\eqn{\sum_j \max\{0,\, b_{j,t} - (F_t^\top \mathbf 1)_j\}}.}
#' }
#'
#' @examples
#' # Ft <- Matrix::Matrix(0, 3, 3, sparse = TRUE)
#' # a_t <- c(10, 5, 0); b_t <- c(6, 6, 3)
#' # external_flux(Ft, a_t, b_t)
#' @export
external_flux <- function(Ft, a_t, b_t) {
  rs <- as.numeric(Matrix::rowSums(Ft))
  cs <- as.numeric(Matrix::colSums(Ft))
  ext_out <- sum(pmax(a_t - rs, 0))
  ext_in  <- sum(pmax(b_t - cs, 0))
  list(ext_out = ext_out, ext_in = ext_in)
}

