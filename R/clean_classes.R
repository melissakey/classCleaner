#' Clean Classes
#' 
#' Test whether each instance in a class actually belongs.
#' @param D A distance matrix containing the pairwise dissimilarity scores between instances
#' @param assignment The assigned group of each instance
#' @param classes The subset of classes on which filtering is performed, or "all" if all classes should be analyzed.
#' @param alpha0 Desired global type I (v1)  or type II (v2) error rate.
#' @param q (v2 only) - the proportion of distances expected to be "close enough" to keep an instance.  Defaults to 0.5.
#' @param labels a vector of labels for each instance.  Must be the same length as D.  If NULL, the algorithm will check for rownames and column names in D.  If none are found, the instances will be labeled with numbers 1:nrow(D).
#' @param exclude_classes names of "mega" classes which should not be included in determining whether or not classCleaner2 is appropriate.  By default, these classes will not be included in the analysis.
#' 
#' @details
#' For each instance in an analyzed class, this function will estimate the probability that it was correctly placed in that class.
#' 
#' @export


classCleaner <- function(D, assignment, classes = "all", alpha0 = 0.05, q = 0.5, labels = NULL, exclude_classes = NULL) {
  
  # Check to make sure distance matrix is a symmetric, non-negative definite matrix and we have an assignment for each entry.
  if(!(is.matrix(D) && isSymmetric(D) && is.numeric(D))) stop("D must be a symmetric matrix with numeric entries")
  if(length(assignment) != nrow(D)) stop("length(asssignment) != nrow(D)")
  if(min(D, na.rm = TRUE) < 0) stop("D should be a distance matrix with entries >= 0.")
  if(length(alpha0) != 1) {
    warning("Multiple values of alpha0 found.  Only the first is used.")
    alpha0 <- alpha0[1]
  }
  
  class_table <- table(assignment)
  
  if(!identical(classes,"all")) {
    Nk <- tryCatch({
      Nk.tmp <- class_table[names(class_table) %in% classes]
    },
      error = function(cond){
        message("An error ocurred in selecting which classes to filter.")
        # message("Here's the original error message:")
        message(cond)
      },
      warning = function(cond){
        message("A problem ocurred in selecting which classes to filter.")
        # message("Here's the original warning message:")
        message(cond)
      })
    ignored_classes <- as.data.frame(class_table[!(names(class_table) %in% classes) | class_table < Nk])
    names(ignored_classes) <- c("class", "no_instances")
    ignored_classes <- within(ignored_classes, {
      v1_no_removed <- NA
      tau <- NA
      t <- NA
      a <- NA
      v1_status <- 'ignored'
    })
  }
  else {
    Nk <- class_table
    ignored_classes <- data.frame()
  }
  
  # class_summary <- as.data.frame(class_table)
  # names(class_summary) <- c('class', 'no_instances')
  # within(class_summary, {
  #   v1_status <- ifelse()
  # })
  
  # remove excluded classes from the list of considered classes
  # save counts for future reporting
  if(!is.null(exclude_classes)) {
    excluded_classes <- as.data.frame(Nk[(names(Nk) %in% exclude_classes)])
    names(excluded_classes) <- c("class", "no_instances")
    excluded_classes <- within(excluded_classes, {
      
      v1_no_removed <- NA
      tau <- NA
      t <- NA
      a <- NA
      v1_status = 'excluded'
    })
    Nk <- Nk[!(names(Nk) %in% exclude_classes)]
  } else   excluded_classes <- data.frame() # empty frame to define variable

  
  # Determine which version(s) of classCleaner are appropriate for the data set
  # v2 is only appropriate if the number of instances in the largest (non-excluded) class is small.
  v2 <- max(Nk / sum(class_table)) < .1
  
  # handle labels
  if(is.null(labels)){
    if(is.null(rownames(D))){
      if(is.null(colnames(D))) labels <- 1:ncol(D)
      else labels <- colnames(D)
    } else labels <- rownames(D)
  }
  
  # If v2 is appropriate, calculate gamma once for each values of N in the data set.
  if(v2) {
    N <- sort(unique(Nk))
    N <- N[N > 1]
    df <- data.frame(
      N = N,
      gamma = vapply(N, function(n) find_p(n - 1, q, alpha = alpha0 / n), 0)
    )
    excluded_classes <- within(excluded_classes, {
      
      v2_no_removed <- NA
      c <- NA
      gamma <- NA
      v2_status <- 'excluded'
    })
    ignored_classes <- within(ignored_classes, {
      
      v2_no_removed <- NA
      c <- NA
      gamma <- NA
      v2_status <- 'ignored'
    })
  }
  
  min_count <- c(20, 2)[v2 + 1]
  result_list <- lapply(names(Nk)[Nk >= min_count], function(k){

    D11 <- D[which(assignment == k), which(assignment == k)]
    D21 <- D[which(assignment != k), which(assignment == k)]
    
    alpha <- alpha0 / Nk[k]
    
    
    psi_t <- psi(D11[lower.tri(D11)], D21)
    
    protein_result <- list(
      Zi = vapply(1:Nk[k], function(i) sum(D11[-i,i] < psi_t["t"]), 0),
      instance = labels[assignment == k],
      index  = which(assignment == k)
    )
    
    
    protein_result <- within(protein_result, {
      if(v2) {
        ind <- df$N == Nk[k]
        tu = vapply(1:Nk[k], function(i) {
          quantile(D[, which(assignment == k)[i]], df$gamma[ind])
        }, 0)
        gamma = df$gamma[ind]
        Ui <- vapply(1:Nk[k], function(i) sum(D11[-i, i] < tu[i]), 0)
        c <- (Nk[k] - 1) * q
        pU <- stats::pbinom(Ui - 1 , Nk[k] - 1, gamma, lower.tail = FALSE)
        pU_bon <- p.adjust(pU, method = 'bonferroni')        
      }
      a <- stats::qbinom(alpha, Nk[k] - 1, psi_t["tau"]) - 1
      tau_hat <- Zi / Nk[k]
      
      p <- stats::pbinom(Zi, Nk[k] - 1, psi_t["tau"])
      p_Bon <- stats::p.adjust(p, method = "bonferroni")
      
      t <- psi_t["t"]
      tau <- psi_t["tau"]
      alpha0 <- alpha0
      Nk <-  as.numeric(Nk[k])
    })
    
    
    v1_summary <- with(protein_result, data.frame(
      class = k,
      v1_status = 'analyzed',
      no_instances = Nk,
      v1_no_removed = sum(p_Bon < 0.05),
      tau = tau,
      t = t,
      a = a
    ))
    rownames(v1_summary) <- NULL
    
    v1_result <- with(protein_result, data.frame(
      class = k,
      instance = instance,
      Zi = Zi,
      p = p,
      p_bonferroni = p_Bon,
      keep = ifelse(p_Bon < 0.05, "remove", "keep")
    ))
    
    if(!v2) return(list(class_summary = v1_summary, v1 = v1_result))
    v2_summary <- with(protein_result, data.frame(
      v2_no_removed = sum(pU_bon >= 0.05),
      c = c,
      gamma = gamma
    ))
  
    v2_result <- with(protein_result, data.frame(
      class = k,
      Ui = Ui,
      t = tu,
      p = pU_bon,
      keep = ifelse(pU_bon < 0.05, "keep", "remove")
    ))
    list(class_summary = rbind(v1_summary, v2_summary), v1 = v1_result, v2 = v2_result, )
  })
  
  result <- list(
    class_summary = rbind(do.call("rbind", lapply(result_list, function(x) {x$class_summary})), ignored_classes, excluded_classes),
    v1_result = do.call("rbind", lapply(result_list, function(x) {x$v1}))
  )
  if(!is.null(result_list[[1]]$v2)) {
    result$v2_result <- do.call("rbind", lapply(result_list, function(x) {x$v2}))
  }
  
  class(result) <- 'classcleaner'
  result
}
