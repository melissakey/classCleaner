find_p <- function(n, q = .5, alpha = 0.05) {
  if(n < 25) return(f1(n, q, alpha))
  else return(f2(n, q, alpha))
}

coefs <- function(n, q, alpha) {
  k <- ceiling(n* q - 1)
  m <- matrix(choose(n, 0:k), ncol = 1)
  
  w <- c(1, -1)
  if(n %% 2) w <- -1 * w
  
  
  mi <- sapply(n - 0:k, function(x) {
    
    if(x %% 2) w <- -1 * w
    
    tmp <- choose(x, 0:x)
    c(rep(0, 1 + n - length(tmp)), tmp) * rep(w, length.out = n + 1)
  })
  
  coefs <- mi %*% m
  coefs[1] <- alpha
  coefs
}

f1 <- function(n, q, alpha, tol = 1e-5) {
  all_ans <- polyroot(coefs(n, q, alpha))
  real_ans <- Re(all_ans[(abs(Im(all_ans) - 0)) < tol])
  real_ans[real_ans > 0 & real_ans < 1]
}


f2 <- function(n , q, alpha) {
  z <- qnorm(1-alpha)
  x <- n + z^2
  
  ret <- (z^2 + 2 * n * q + c(-1, 1) * z *  sqrt(z^2 + 4 * n * q * (1 - q))) / 2 / (z^2 + n)
  
  min(ret)
}
