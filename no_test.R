no_test <- function(x,y,qi=c(0.05,0.5,0.95), nboot=1999){
  # Compares two independent groups through different quantiles based on 
  # Mahalanobis distance and a percentile bootstrap approach.
  # The NO quantile estimator is used.
  # requires noquantile() function
  dif <- noquantile(x,qi)-noquantile(y,qi)
  bvec <- matrix(NA,nrow=nboot,ncol=length(qi))
  n <- sapply(list(x,y), length)
  for(i in 1:nboot){
    xi <- sample(x,size=n[1],replace=T)
    yi <- sample(y,size=n[2],replace=T)
    bvec[i,] <- (noquantile(xi,qi)-noquantile(yi,qi))
  }
  vecz <- rep(0,length(qi))
  smat <- var(bvec)
  bvec <- rbind(bvec,vecz)
  db <- sqrt(mahalanobis(bvec,dif,smat))
  sig.level <- sum(db[nboot+1] < db[1:nboot])/nboot
  sig.level
}