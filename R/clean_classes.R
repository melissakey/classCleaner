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

  # handle labels
  if(is.null(labels)){
    if(is.null(rownames(D))){
      if(is.null(colnames(D))) labels <- 1:ncol(D)
      else labels <- colnames(D)
    } else labels <- rownames(D)
  }
    
  class_table <- table(assignment)
  
  if(identical(classes, "all")) {
    classes <- names(class_table)
  }
  else {
    unavailable_classes <- setdiff(classes, names(class_table))
    if(length(unavailable_classes) > 0) stop({cat("Some reported classes not found.\n"); print(unavailable_classes )})
  }
  
  class_summary <- as.data.frame(class_table, row.names = names(class_table))
  names(class_summary) <- c('class', 'no_instances')

  class_summary <- within(class_summary, {
    v1_status <- lest::case_when(
      class %in% exclude_classes ~ 'excluded',
      class %in% classes & no_instances >= 20 ~ 'analyzed',
      TRUE ~ 'ignored'
    )
    v1_no_removed <- NA
    tau <- NA
    t_star <- NA
    a <- NA   
  })
  
  # Determine which version(s) of classCleaner are appropriate for the data set
  # v2 is only appropriate if the number of instances in the largest (non-excluded) class is small.
  v2 <- with(class_summary, max(no_instances[v1_status != 'excluded'] / sum(no_instances)) < .1)

 
  # If v2 is appropriate, calculate gamma once for each values of N in the data set.
  if(v2) {
    class_summary <- within(class_summary, {
      v2_status <- lest::case_when(
        class %in% exclude_classes ~ 'excluded',
        class %in% classes & no_instances > 1 ~ 'analyzed',
        TRUE ~ 'ignored'
      )
      v2_no_removed <- NA
      c <- NA
      gamma <- NA
    })
    
    N <- with(class_summary, sort(unique(no_instances[v2_status == 'analyzed'])))
    df <- data.frame(
      N = N,
      gamma = vapply(N, function(n) find_p(n - 1, q, alpha = alpha0 / n), 0)
    )
    
  }
  analyze_list <- with(class_summary, levels(class[v1_status == "analyzed" | (v2 && v2_status == "analyzed")]))
  
  
  result_list <- lapply(analyze_list, function(k){
    Nk <- with(class_summary, no_instances[class == k])
    alpha <- alpha0 / Nk
    
    D11 <- D[which(assignment == k), which(assignment == k)]    
    
    # indexing
    instance = labels[assignment == k]
    index  = which(assignment == k)
    
    # v1
    if(class_summary[k, "v1_status"] == 'analyzed') {
      D21 <- D[which(assignment != k), which(assignment == k)]
      
      psi_t <- psi(D11[lower.tri(D11)], D21)
      Zi = vapply(1:Nk, function(i) sum(D11[-i,i] < psi_t["t"]), 0)
      p1 <- stats::pbinom(Zi, Nk - 1, psi_t["tau"])
      p1_bon <- stats::p.adjust(p1, method = "bonferroni")
      
      
      class_summary[k, 'a'] <<- stats::qbinom(alpha, Nk - 1, psi_t["tau"]) - 1
      class_summary[k, 't_star'] <<- psi_t['t']
      class_summary[k, 'tau'] <<- psi_t['tau']
      class_summary[k, 'v1_no_removed'] <<- sum(p1_bon < 0.05)
      
      v1_result <- data.frame(
        class = k,
        instance = instance,
        Zi = Zi,
        p = p1,
        p_bonferroni = p1_bon,
        keep = ifelse(p1_bon < 0.05, "remove", "keep")
      )
    } else v1_result <- data.frame()

    # v2
    if(v2 && class_summary[k, "v2_status"] == 'analyzed') {
      ind <- df$N == Nk
      
      t_star_star = vapply(index, function(i) {quantile(D[, i], df$gamma[ind])}, 0)
      Ui <- vapply(1:Nk, function(i) sum(D11[-i, i] < t_star_star[i]), 0)
      
      p2 <- stats::pbinom(Ui - 1 , Nk - 1,  df$gamma[ind], lower.tail = FALSE)
      p2_bon <- p.adjust(p2, method = 'bonferroni')
      
      class_summary[k, 'c'] <<- (Nk - 1) * q
      class_summary[k, 'gamma'] <<-  df$gamma[ind]
      class_summary[k, 'v2_no_removed'] <<- sum(p2_bon >= 0.05)
      
      v2_result <- data.frame(
        class = k,
        Ui = Ui,
        t_star2 = t_star_star,
        p = p2_bon,
        keep = ifelse(p2_bon < 0.05, "keep", "remove")
      )
      
    } else v2_result <- data.frame()
    
    list(v1 = v1_result, v2 = v2_result)
  })
  

  result <- list(
    class_summary = class_summary,
    v1_result = do.call("rbind", lapply(result_list, function(x) {x$v1})),
    v2_result = do.call("rbind", lapply(result_list, function(x) {x$v2}))
  )

  class(result) <- 'classcleaner'
  result
}
