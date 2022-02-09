#' Simulate clustered data using the normal distribution.
#'
#' This function generates a matrix with \eqn{n} indepenent rows and \eqn{N} columns,
#' where the columns are clustered into \eqn{K} classes with \eqn{N_k}{N[k]} instances in class \eqn{k}
#' @param n The total number of observations per instance.
#' @param Nk A vector of length \eqn{K} giving the number of  within each class
#' @param s A vector of standard deviations. See details.
#' @param rho A matrix of correlation coefficients.  See details.
#' @param tau The within-group variance.  Only used when method = "by_class".
#' @param method Either "by_instance" or "by_class".
#' @details
#' This function generates a matrix with a block - correlation stucture across columns and independent rows.
#' When method = 'by_instance', the values of s and rho are taken to be instance-level properties of the data.
#' That is, s is a vector of length K such that the ith entry is the standard deviations of observations within class k
#' and rho is a  \eqn{K \times K}{K * K} symmetric matrix such that entry (i,j) gives the correlation between an instance is class i and an instance in class j.  Correspondingly, entry (i,i) gives the correlation between two (different) instances in class i.
#' In contrast, when method = 'by_class', the values of s and rho are taken to be class-level properties.
#' The variance from the 'by_instance' characterization is broken down into a class-level variance (s^2), which gives the variability of the "true" pattern of the class over the observations"  and an instance-level variance (tau^2).which gives the variabilty of the the observed instances from the true pattern. The correlation is now in terms of the classes: entry (i,j) gives the correlation between the "true" pattern of class i and the "true" pattern of class j.
#' @export
#' @examples
#' rho.mat <- matrix(c(.5,.2,.2,.3),nrow = 2,ncol = 2)
#' X <- simulate_by_instance(Nk = c(50, 100),rho = rho.mat,n = 150)

simulate_clustered_data <- function(
  n = 100,        # total number of observations per instance
  Nk = c(40,200), # number of instances in each group
  s = c(1,1),
  rho = matrix(c(.6,.1,.1,.25),nrow = 2,ncol = 2),
  tau = 1,
  method = c("by-class", "by-instance")
) {
  method = match.arg(method)

  # error checking
  if (length(s) == 1) s <- rep(s, length(Nk))
  if(!is.matrix(rho)) rho <- as.matrix(rho)
  
  if (nrow(rho) != ncol(rho)) stop("rho must be a square matrix.")
  
  
  if (length(Nk) != length(s)) stop("s must be of length 1 or have the same length as Nk")
  if (length(tau) > 1) stop ("support for tau > 1 is not yet implemented.")
  if (!all.equal(rho[lower.tri(rho)], rho[upper.tri(rho)]))
    warning("rho is assumed to be a symmetric matrix.  Only the lower triangle is used.")

  # compute it
  if(method == "by-instance"){
    if (length(Nk) != nrow(rho)) stop("Nk must have the same length as ncol(rho)/nrow(rho)")
    X <- sim_by_instance(n, Nk, s, rho)
  }
  else{
    if (length(Nk) == nrow(rho)){
      X <- sim_by_class(n, Nk[Nk > 0], s[Nk > 0], tau, rho[Nk > 0, Nk > 0])
    } else if (nrow(rho) == 1) {
      X <- sim_by_class(n, Nk[Nk > 0], s[Nk > 0], tau, as.matrix(rho))
    } else       stop("Nk must have the same length as ncol(rho)/nrow(rho) or rho must be of length 1.")

    
  }

  # add identifiers
  rownames(X) <- paste("n",1:n,sep = ".")
  colnames(X) <- paste0("N",rep(1:length(Nk),Nk),".",
    unlist(sapply(Nk, function(x) {if(x > 0) 1:x else integer(0)})))

  X
}

