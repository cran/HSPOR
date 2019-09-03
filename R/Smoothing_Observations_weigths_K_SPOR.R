
.Smoothing_Observations_weigths_K_SPOR <- function(X,deg,timeT){

  nbT <- length(timeT)
  obs_weights <- matrix(0,length(X),nbT+1)
  m <- 1

  while(m <= (nbT)){

    e1 = X[(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1))]
    e2 = X[(which(X > timeT[m])[1] + (deg))]

    if(m == 1){

      if(m == length(timeT)){
        pos1 <- min(1,(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)):(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)
        pos2 <-(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)):(which(X > timeT[m])[1] + (deg))
        pos3 <- (which(X > timeT[m])[1] + (deg) + 1):length(X)
        pos = c(pos1,pos2,pos3)

        p1 = ((-1/(e2-e1)^2)*(X[pos2] - e1)^2)+1
        if(pos1[1] == 0){
          if(pos3[1] == 0){
            pos = pos2
            vz1 <- p1
            vz2 <- 1-p1
          }else{
            pos <- c(pos2,pos3)
            vz1 <-c(p1,rep(0,length(pos3)))
            vz2 <- c(1-p1,rep(1,length(pos3)))
          }
        }else{
          if(pos3[1] == 0){
            pos <- c(pos1,pos2)
            vz1 <- c(rep(1,length(pos1)),p1)
            vz2 <- c(rep(0,length(pos1)),1-p1)
          }else{
            pos <- c(pos1,pos2,pos3)
            vz1 <-c(rep(1,length(pos1)),p1,rep(0,length(pos3)))
            vz2 <- c(rep(0,length(pos1)),1-p1,rep(1,length(pos3)))
          }
        }

        obs_weights[pos,m] <- vz1
        obs_weights[pos,m+1] <-vz2
      } else{
        pos1 <- min(1,(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)):(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)
        pos2 <-(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)):(which(X > timeT[m])[1] + (deg))

        p1 = ((-1/(e2-e1)^2)*(X[pos2] - e1)^2)+1
        if(pos1[1] == 0){
          pos = pos2
          vz1 <- p1
          vz2 <- 1-p1
        }else{
          pos = c(pos1,pos2)
          vz1 <- c(rep(1,length(pos1)),p1)
          vz2 <- c(rep(0,length(pos1)),1-p1)
        }

        obs_weights[pos,m] <- vz1
        obs_weights[pos,m+1] <-vz2
      }
    }

    if(m > 1 & m < (nbT+1)){

      if(m == length(timeT)){
        pos1 <- (which(X > timeT[m-1])[1] + (deg) + 1):(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)
        pos2 <- (which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)):(which(X > timeT[m])[1] + (deg))
        pos3 <- (which(X > timeT[m])[1] + (deg) + 1):length(X)

        p2 = ((-1/(e2-e1)^2)*(X[pos2] - e1)^2)+1
        if(pos1[1] == 0){
          if(pos3[1] == 0){
            pos = pos2
            vz3 <- p2
            vz4 <- 1-p2
          }else{
            pos = c(pos2,pos3)
            vz3 <-c(p1,rep(0,length(pos3)))
            vz4 <- c(1-p2,rep(1,length(pos3)))
          }
        }else{
          if(pos3[1] == 0){
            pos = c(pos1,pos2)
            vz3 <- c(rep(1,length(pos1)),p2)
            vz4 <- c(rep(0,length(pos1)),1-p2)
          }else{
            pos = c(pos1,pos2,pos3)
            vz3 <-c(rep(1,length(pos1)),p2,rep(0,length(pos3)))
            vz4 <- c(rep(0,length(pos1)),1-p2,rep(1,length(pos3)))
          }
        }

        obs_weights[pos,m] <- vz3
        obs_weights[pos,m+1] <-vz4
      }
      else{
        pos1 <- (which(X > timeT[m-1])[1] + (deg) + 1):(which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)-1)
        pos2 <- (which(X <= timeT[m])[length(which(X <= timeT[m]))] - (deg+1)):(which(X > timeT[m])[1] + (deg))

        p2 = ((-1/(e2-e1)^2)*(X[pos2] - e1)^2)+1

        if(pos1[1] == 0){
          pos = pos2
          vz3 <- p2
          vz4 <- 1-p2
        }else{
          pos = c(pos1,pos2)
          vz3 <- c(rep(1,length(pos1)),p2)
          vz4 <- c(rep(0,length(pos1)),1-p2)
        }

        obs_weights[pos,m] <- vz3
        obs_weights[pos,m+1] <-vz4
      }
    }
    m <- m+1
  }

  return(obs_weights)
}





