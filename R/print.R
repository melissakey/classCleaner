#' @export
print.classcleaner <- function(cc_obj, sort = TRUE, ...) {

  K <- nrow(cc_obj$class_summary)

  # check to see if n is defined
  args <- list(...)
  if(is.null(args$n)) n <- if(K > 20) 10
  
  
  v2 <- !is.null(cc_obj$v2_result)
  
  sum_obj <- cc_obj$class_summary
  rownames(sum_obj) <- sum_obj$class
  sum_obj <- within(sum_obj, {
    v1_no_retained <- no_instances - v1_no_removed
    pcnt <- round(v1_no_retained / no_instances * 100, digits = 1)
  })
  
  sum_obj <- sum_obj[, c('no_instances', 'v1_no_retained', 'pcnt')]
  names(sum_obj) <- c("# Instances", "# Retained", "%")

  cat(paste("classCleaner was applied to", nrow(sum_obj), "classes\n" ))
  cat(paste("A total of", sum(sum_obj$`# Retained`), "instances were retained out of", sum(sum_obj$`# Instances`),"\n\n"))
  if(sort)  {
    sum_obj <- sum_obj[order(sum_obj$`# Retained`, decreasing = TRUE),]
  }
  print(sum_obj, n, ...)
    
  if(v2) {
    sum_obj <- cc_obj$class_summary
    rownames(sum_obj) <- sum_obj$class
    sum_obj <- within(sum_obj, {
      v2_no_retained <- no_instances - v2_no_removed
      pcnt <- round(v2_no_retained / no_instances * 100, digits = 1)
    })
    
    sum_obj <- sum_obj[, c('no_instances', 'v2_no_retained', 'pcnt')]
    names(sum_obj) <- c("# Instances", "# Retained", "%")
    
    cat(paste("classCleaner2 was applied to", nrow(sum_obj), "classes\n" ))
    cat(paste("A total of", sum(sum_obj$`# Retained`), "instances were retained out of", sum(sum_obj$`# Instances`),"\n\n"))
    if(sort)  {
      sum_obj <- sum_obj[order(sum_obj$`# Retained`, decreasing = TRUE),]
    }
    print(sum_obj, n, ...)
  }
}
