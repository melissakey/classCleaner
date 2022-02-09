#' @export
print.classcleaner <- function(cc_obj, sort = TRUE, ...) {

  K <- nrow(cc_obj$class_summary)
  
  # check to see if n is defined
  args <- list(...)
  if(is.null(args$n)) {
    n <- if(K > 20) 10 else K
  } else n <- args$n
  
  if(nrow(cc_result$v1_result) > 0) {
    sum_obj <- cc_obj$class_summary
    rownames(sum_obj) <- sum_obj$class
    sum_obj <- within(sum_obj, {
      v1_no_retained <- no_instances - v1_no_removed
      pcnt <- round(v1_no_retained / no_instances * 100, digits = 1)
    })
    
    sum_obj <- sum_obj[, c('no_instances', 'v1_no_retained', 'pcnt')]
    names(sum_obj) <- c("# Instances", "# Retained", "%")
    
    Ka <- sum(cc_obj$class_summary$v1_status == 'analyzed')
    
    cat(paste("classCleaner was applied to", Ka, "of", K, "classes\n"))
    cat(paste("A total of", sum(sum_obj$`# Retained`, na.rm = TRUE), "instances were retained out of", sum(sum_obj$`# Instances`),"\n"))
    if(sum(cc_obj$class_summary$no_instances < 20) > 0)
      cat(paste("A total of", sum(cc_obj$class_summary$no_instances < 20), "instances were ignored due to small sample size.\n\n"))
    if(sort)  {
      sum_obj <- sum_obj[order(sum_obj$`# Retained`, decreasing = TRUE),]
    }
    print(head(sum_obj, n), ...)
  }
  if((nrow(cc_result$v1_result) > 0) & nrow(cc_result$v2_result) > 0) cat("\n\n")
  if(nrow(cc_result$v2_result) > 0) {
    sum_obj <- cc_obj$class_summary
    rownames(sum_obj) <- sum_obj$class
    sum_obj <- within(sum_obj, {
      v2_no_retained <- no_instances - v2_no_removed
      pcnt <- round(v2_no_retained / no_instances * 100, digits = 1)
    })
    
    sum_obj <- sum_obj[, c('no_instances', 'v2_no_retained', 'pcnt')]
    names(sum_obj) <- c("# Instances", "# Retained", "%")
    
    Ka <- sum(cc_obj$class_summary$v2_status == 'analyzed')
    
    cat(paste("classCleaner2 was applied to", Ka, "of", K, "classes\n"))
    cat(paste("A total of", sum(sum_obj$`# Retained`, na.rm = TRUE), "instances were retained out of", sum(sum_obj$`# Instances`),"\n\n"))
    # cat(paste("A total of", sum(cc_obj$class_summary$no_instances < 2), "instances were ignored due to small sample size.\n\n"))
    if(sort)  {
      sum_obj <- sum_obj[order(sum_obj$`# Retained`, decreasing = TRUE),]
    }
    print(head(sum_obj, n), ...)
  }
}

#' @export
predict.classcleaner <- function(cc_obj, method = c('all-available', 'v1', 'v2')) {
  method <- match.arg(method)
  K1 <- nrow(cc_obj$v1_result)
  K2 <- nrow(cc_obj$v2_result)
  
  if(method == 'all-available' & K1 > 0 & K2 > 0) {
    df <- merge(cc_obj$v1_result, cc_obj$v2_result, all = TRUE, by = c("class", "instance"), suffixes = c(".v1", '.v2'))[, c('class', 'instance', 'keep.v1', 'keep.v2')]
  } else if((method == 'v1' & K1 > 0) || (method == 'all-available' & K1 > 0)) {
    df <- cc_obj$v1_result[, c('class', 'instance', 'keep')]
  } else if((method == 'v2' & K2 > 0) || (method == 'all-available' & K2 > 0)) {
    df <- cc_obj$v2_result[, c('class', 'instance', 'keep')]
  } else {
    stop("Could not determine which results to predict")
  }
  df
}