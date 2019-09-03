
.Minus_Complete_log_likelihood <- function(X,Y,deg,timeT,mat_param,sigma2,obs_weights){

  r <- 1
  wp <- 0
  nbT <- length(timeT)
  while(r <= (nbT+1)){
    s <- 0
    for(i in 1:(deg+1)){
      s <- s + mat_param[r,i] * X^(deg+1-i)
    }

    wp <- wp + obs_weights[,r]*((1/sqrt(2*pi*sigma2[r])) * exp( - ((Y - s)^2)/(2*sigma2[r])))

    r <- r+1
  }

  Complete_log_likelihood <- sum(log(wp))

  return(-Complete_log_likelihood)
}
