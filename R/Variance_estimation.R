
.Variance_estimation <- function(X,Y,deg,timeT,mat_param,regimes_prob){

  nbT <- length(timeT)

  r <- 1
  sigma2 <- matrix(0,nrow=1,ncol=nbT+1)
  while(r <= (nbT+1)){
    s <- 0
    for(i in 1:(deg+1)){
      s <- s + mat_param[r,i] * X^(deg+1-i)
    }
    sigma2[,r] <- sum(((Y - s)^2) * regimes_prob[,r])/sum(regimes_prob[,r])
    r <- r+1
  }

  return(sigma2)
}


