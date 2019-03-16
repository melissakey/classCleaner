#' Generate data
#' A convience wrapper around simulate_clustered_data for generating clustered data sets
#' 
#' @param Nk The number of instances in each class (can be a table)
#' @param p_mismatch The expected proprotion of instances in any given class which are misclassified
#' @param p_garbage The expected proportion of instances in any given class which are actually quantitative noise.
#' @param n The number of samples
#' @param s The within-class standard deviation
#' @param tau The between-class standard deviation
#' @param rho The average correlation between different classes (must be greater than -1 / K, where K is the number of classes)
#' 
#' @export

generate_data <- function(Nk, p_mismatch, p_garbage, K_proteome, n, s = 1, tau = 1, rho = 0.0) {
  
  ######-------------------------------------------------------------------------
  #
  # Set internal parameters
  #
  ######-------------------------------------------------------------------------
  
  K <- length(Nk)
  N <- sum(Nk)
  
  correct_assignment <- rep(1:K, Nk)
  U1 <- runif(N)
  U2 <- rbeta(N, 2 / p_mismatch - 2, 2)
  
  alternate_assigned_protein <- sample(K_proteome, size = N, replace = TRUE, prob = c(Nk * 10, rep(1, K_proteome - K)))
  actual_assignment <- ifelse(U2 > U1, correct_assignment, alternate_assigned_protein)
  confidence <- pmax(U1, U2) * 100
  
  ######-------------------------------------------------------------------------
  #
  # Generate the "real" abundance data
  #
  ######-------------------------------------------------------------------------
  
  X <- tryCatch(
    classCleaner::simulate_clustered_data(n = n, Nk = Nk, s = rep(s, length(Nk)), tau = tau, rho = 0.0, method = "by-class"),
    error = function(cond) {
      # print(paste0("Failed to generate data for parameter set ",set_num,", simulation ",b))
      return(NULL)
    }
  )
  
  if(is.null(X)) return(NULL)
  
  ######-------------------------------------------------------------------------
  #
  # Replace with garbage if needed
  #
  ######-------------------------------------------------------------------------
  
  N_garbage <- floor(N * p_garbage)
  
  if (N_garbage > 0) {
    garbage_ind <- sample(N, N_garbage)                     # pick values to corrupt
    X_garbage <- classCleaner::simulate_clustered_data(n = n, Nk = rep(1,N_garbage),s = 1, tau = 1, rho = 0, method = "by-class")
    X[, garbage_ind] <- X_garbage - .5
    correct_assignment[garbage_ind] <- 0
    colnames(X)[garbage_ind]  <- paste0(stringr::str_replace_all(colnames(X)[garbage_ind], "(?<=N)\\d+.*",c("0.")), 1:length(garbage_ind))
  } else {
    garbage_ind = NULL
  }
  
  
  list(
    actual_assignment = actual_assignment,
    correct_assignment = correct_assignment,
    garbage_indices = garbage_ind,
    confidence = confidence,
    X = X
  )
}
