#' classCleaner: A package for cleaning outliers when data is grouped.
#' 
#' @docType package
#' @name classCleaner
#' 
#' @examples 
#' set.seed(23)
#' 
#' X <- simulate_clustered_data(
#'   n = 200,
#'   Nk = rep(50, 100),
#'   s = rep(1, 100),
#'   rho = .2,
#'   tau = 1,
#'   method = "by-class"
#' )
#' # true assignment
#' a <- rep(1:100, each = 50)
#' 
#' # corrupted assignment
#' b <- sample(100, 50 * 100, replace = TRUE)
#' 
#' # corrupt 10% of samples
#' a.corrupt <- ifelse(runif(50 * 1000) < 0.1, b, a)
#' 
#' D <- 1 - cor(X)
#' result <- identify_outliers(a.corrupt, D, 1000, colnames(D))
#' 
#' @useDynLib classCleaner, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @improt RcppProgress
NULL