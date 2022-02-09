psi_old <-function(D1, D2){
  F1<-stats::ecdf(D1)
  F2<-stats::ecdf(D2)
  N<-4000
  tt<-seq(min(D1, na.rm = TRUE),max(D2, na.rm = TRUE), length=N)
  
  delta <- abs(1 - F2(tt) - F1(tt))
  tc <-  stats::median(tt[which(delta <= delta[which.min(delta)] + 10e-8)])
  c(
    t = tc,
    tau = F1(tc)
  )
}
