
.Variance_estimation_K_SPOR_DynProg <- function(X,Y,deg,mat_param){

  s <- 0
  for(i in 1:(deg+1)){
    s <- s + mat_param[1,i] * X^(deg+1-i)
  }
  sigma2 <- sum(((Y - s)^2))/length(X)

  return(sigma2)
}


