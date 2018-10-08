#' Clean Classes
#' 
#' Test whether each instance in a class actually belongs.
#' @param D A distance matrix containing the pairwise dissimilarity scores between instances
#' @param assignment The assigned group of each instance
#' @param classes The subset of classes on which filtering is performed, or "all" if all classes should be analyzed.
#' @param labels identifier for each instance.  If NULL, the row/column names or indices of D are used
#' 
#' @details
#' For each instance in an analyzed class, this function will estimate the probability that it was correctly placed in that class.
#' 
#' @export


clean_classes <- function(D, assignment, classes = 'all', tau = c(0.01, 0.05, 0.1), labels = NULL, display_progress = FALSE) {
  
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
      Nk.tmp <- Nk[classes]
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
  
  # handle labels
  if(is.null(labels)){
    if(is.null(rownames(D))){
      if(is.null(colnames(D))) labels <- 1:ncol(D)
      else labels <- colnames(D)
    } else labels <- rownames(D)
  }
  
  
  result <- lapply(names(Nk)[Nk > 1], {
    D11 <- D[which(assignment == k), which(assignment == k)]
    D21 <- D[which(assignment == k), which(assignment != k)]
    
    
    
    tc <- quantile(D21, tau)
    
    Zi <- as.data.frame(t(sapply(1:Nk[k], function(i) colSums(outer(D11[-i, i], tc, "<")))))
    
    Zi <- reshape(Zi,
      direction = 'long',
      idvar = 'instance',
      ids = labels[assignment == k],
      timevar = 'cutoff',
      times = names(Zi),
      v.names = "Z",
      varying = list(names(Zi)),
      # new.row.names = 1:prod(dim(Zi))
    )
    
    Zi <- within(Zi, {
      tau <- as.numeric(sub("%","",cutoff)) / 100
      Nk <- as.numeric(Nk[k])
      
      p <- 1 - pbinom(Z, Nk - 1, tau)
      
      k <-  k
      index <- match(instance, labels)
      
      # if(!is.null(colnames(D11))) label <- colnames(D11)
      # else label <- paste0("N",k,".i",1:Nk[k])
      
      
      p_BH <- p.adjust(p, method = "BH")
      p_BY <- p.adjust(p, method = 'BY')
    })
  })
  result <- do.call("rbind", result)
  result[c("k", "index", "label", 'tau', "Zi", "p", "p_BH", "p_BY", "Nk")]
}