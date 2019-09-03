
.Observations_weigths <- function(X,timeT){

  nbT <- length(timeT)
  pZ <- matrix(0,length(X),nbT+1)
  m <- 1
  while(m <= (nbT +1)){
    if(m == 1){
      pos <- which(X <= timeT[m])
    }
    if(m > 1 & m < (nbT+1)){
      pos <- which(X > timeT[m-1] & X <= timeT[m])
    }
    if(m == (nbT+1)){
      pos <- which(X > timeT[m-1])
    }
    vZ <- rep(1,length(pos))
    pZ[pos,m] <- vZ
    m <- m+1
  }
  return(pZ)
}
