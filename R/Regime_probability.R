
.Regime_probability <- function(X,Y,deg,timeT,mat_param,sigma2,obs_weights){

  nbT <- length(timeT)

  r <- 1
  regimes_prob <- matrix(0,length(X),nbT+1)
  w <- matrix(0,length(X),nbT+1)
  f <- 0

  while(r <= (nbT+1)){
    s <- 0

    for(i in 1:(deg+1)){
      s <- s + mat_param[r,i] * X^(deg+1-i)
    }

    w[,r] <- (1/sqrt(2*pi*sigma2[r])) * exp( - ((Y - s)^2)/(2*sigma2[r]))

    f <- f + (obs_weights[,r] * w[,r])

    r <- r+1
  }

  for(l in 1:(nbT+1)){
    regimes_prob[,l] <- (obs_weights[,l]*w[,l])/f
    if(length(which(is.na(regimes_prob[,l]))) > 0){
      regimes_prob[which(is.na(regimes_prob[,l])),l] <- round(obs_weights[which(is.na(regimes_prob[,l])),l])
    }
  }


  return(regimes_prob)
}



