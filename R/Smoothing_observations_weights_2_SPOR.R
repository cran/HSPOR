
.Smoothing_observations_weigths_2_SPOR <- function(X,deg,timeT){

  nbT <- length(timeT)
  obs_weights <- matrix(0,length(X),nbT+1)

  e1 = X[(which(X <= timeT)[length(which(X <= timeT))] - (deg+1))]
  e2 = X[(which(X > timeT)[1] + (deg))]

  pos1 <- min(1,(which(X <= timeT)[length(which(X <= timeT))] - (deg+1)-1)):(which(X <= timeT)[length(which(X <= timeT))] - (deg+1)-1)
  pos2 <-(which(X <= timeT)[length(which(X <= timeT))] - (deg+1)):(which(X > timeT)[1] + (deg))
  pos3 <- (which(X > timeT)[1] + (deg) + 1):length(X)

  p1 = ((-1/(e2-e1)^2)*(X[pos2] - e1)^2)+1
  if(pos1[1] == 0){
    if(pos3[1] == 0){
      vz1 <- p1
      vz2 <- 1-p1
    }else{
      vz1 <-c(p1,rep(0,length(pos3)))
      vz2 <- c(1-p1,rep(1,length(pos3)))
    }
  }else{
    if(pos3[1] == 0){
      vz1 <- c(rep(1,length(pos1)),p1)
      vz2 <- c(rep(0,length(pos1)),1-p1)
    }else{

      vz1 <-c(rep(1,length(pos1)),p1,rep(0,length(pos3)))
      vz2 <- c(rep(0,length(pos1)),1-p1,rep(1,length(pos3)))
    }
  }
  obs_weights[,1] <- vz1
  obs_weights[,2] <- vz2

  return(obs_weights)
}
