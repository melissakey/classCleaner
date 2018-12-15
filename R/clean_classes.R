#' Clean Classes
#' 
#' Test whether each instance in a class actually belongs.
#' @param D A distance matrix containing the pairwise dissimilarity scores between instances
#' @param assignment The assigned group of each instance
#' @param classes The subset of classes on which filtering is performed, or "all" if all classes should be analyzed.
#' @param labels identifier for each instance.  If NULL, the row/column names or indices of D are used
#' @param tau cutoff for F distribution (using F-distribution only)
#' @param alpha0 desired overall, FWER-style type I error rate
#' @param beta0 desired overall FNR-style type II error rate 
#' @param labels a vector of labels for each instance.  Must be the same length as D.  If NULL, the algorithm will check for rownames and column names in D.  If none are found, the instances will be labeled with numbers 1:nrow(D).
#' @param display_progress future progress bar (not used yet)
#' 
#' @details
#' For each instance in an analyzed class, this function will estimate the probability that it was correctly placed in that class.
#' 
#' @export


clean_classes <- function(D, assignment, classes = 'all', alpha0 = 0.05, beta0 = 0.05, labels = NULL, display_progress = FALSE) {
  
  # Check to make sure distance matrix is a symmetric, non-negative definite matrix and we have an assignment for each entry.
  if(!(is.matrix(D) && isSymmetric(D) && is.numeric(D))) stop('D must be a symmetric matrix with numeric entries')
  if(length(assignment) != nrow(D)) stop('length(asssignment) != nrow(D)')
  if(min(D) < 0) stop("D should be a distance matrix with entries >= 0.")
  
  class_table <- table(assignment)
  
  # progressbar
  if(display_progress)
    pb <- utils::txtProgressBar(min = 0, max = length(class_table), style = 3)
  
  if(!identical(classes,"all")) {
    Nk <- tryCatch({
      Nk.tmp <- class_table[classes]
    },
      error = function(cond){
        message("An error ocurred in selecting which classes to filter.")
        message("Here's the original error message:")
        message(cond)
      },
      warning = function(cond){
        message("A problem ocurred in selecting which classes to filter.")
        message("Here's the original warning message:")
        message(cond)
      })
  }
  else {Nk <- class_table}

  # handle labels
  if(is.null(labels)){
    if(is.null(rownames(D))){
      if(is.null(colnames(D))) labels <- 1:ncol(D)
      else labels <- colnames(D)
    } else labels <- rownames(D)
  }
  
  result <- lapply(names(Nk)[Nk > 1], function(k){
   
    D11 <- D[which(assignment == k), which(assignment == k)]
    D21 <- D[which(assignment == k), which(assignment != k)]
    
    alpha <- alpha0 / Nk[k]


    psi_t <- psi(D11[lower.tri(D11)], D21)
    # m_tc <- c(quantile(D21, tau), psi_t['t'])
    # names(m_tc)[length(m_tc)] <- paste0(psi_t['tau'] * 100, "%")
    # pp <- sapply(1:Nk[k], function(i) colMeans(outer(D21[i,], m_tc, "<")))
    
    Zi_psi <- data.frame(
      Zi = vapply(1:Nk[k], function(i) sum(D11[-i,i] < psi_t['t']), 0),
      instance = labels[assignment == k],
      index  = which(assignment == k)
    )
    Zi_psi <- within(Zi_psi[order(Zi_psi$Zi),], {
      ca <- stats::qbinom(alpha, Nk[k] - 1, psi_t['tau'])
      cb <- stats::qbinom(alpha, Nk[k] - 1, 1 - psi_t['tau'], lower.tail = FALSE) 
      q <- 1 - stats::pbinom(Zi, Nk[k] - 1, 1 - psi_t['tau'])
      q_BH <- stats::p.adjust(q, method = 'BH')
      q_BY <- stats::p.adjust(q, method = 'BY')
      q_Bon <- stats::p.adjust(q, method = 'bonferroni')

      p <- stats::pbinom(Zi, Nk[k] - 1, psi_t['tau'])
      p_BH <- stats::p.adjust(p, method = 'BH')
      p_BY <- stats::p.adjust(p, method = 'BY')
      p_Bon <- stats::p.adjust(p, method = 'bonferroni')

      t <- psi_t['t']
      tau <- psi_t['tau']
      alpha0 <- alpha0
      Nk <-  as.numeric(Nk[k])
      k <- factor(k)
    })
  })
  result <- do.call("rbind", result)
}