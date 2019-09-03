
.Jacobian_Matrix_K_SPOR <- function(X,Y,deg,regimes_prob,timeT,sigma2,constraint){

  nbT = length(timeT)

  d = deg
  if(constraint == 0){
    t_mat <- (nbT+1)*(d+1)
  }else if(constraint == 1){
    t_mat <- nbT*(d+2) + (d+1) #nbT constraints parameters
  }else{
    t_mat <- (nbT+1)*(d+1) + (nbT*2) #nbT*2 constraints parameters
  }

  if(constraint == 1){ #if continuity constraint
    m1 <- matrix(0,nbT,t_mat)
    m2 <- matrix(0,(nbT+1)*(d+1),nbT)
    i=1
    while(i <= nbT){
      vF1 <- c()
      vF2 <- c()
      deg <- d

      while(deg >= 0){
        v1 <- rep(0,nbT+1)
        v2 <- rep(0,nbT+1)
        v1[c(i,i+1)] <- c(timeT[i]^deg,-timeT[i]^deg)
        v2[c(i,i+1)] <- c((timeT[i]^deg)*sigma2[i],-(timeT[i]^deg)*sigma2[i+1])
        vF1 <- c(vF1,v1)
        vF2 <- c(vF2,v2)
        deg=deg-1
      }
      m1[i,] <- c(vF1,rep(0,nbT))
      m2[,i] <- vF2
      i <- i+1
    }
  }

  if(constraint == 2){ #if differentiability constraint
    m11 <- matrix(0,nbT,t_mat)
    m12 <- matrix(0,nbT,t_mat)
    m21 <- matrix(0,(nbT+1)*(d+1),nbT)
    m22 <- matrix(0,(nbT+1)*(d+1),nbT)
    i<-1
    while(i <= nbT){
      vc1 <- c()
      vc2 <- c()
      vd1 <- c()
      vd2 <- c()
      deg <- d
      while(deg >= 0){

        #continuity constraint
        v1 <- rep(0,nbT+1)
        v2 <- rep(0,nbT+1)
        v1[c(i,i+1)] <- c(timeT[i]^deg,-timeT[i]^deg)
        v2[c(i,i+1)] <- c((timeT[i]^deg)*sigma2[i],-(timeT[i]^deg)*sigma2[i+1])
        vc1 <- c(vc1,v1)
        vc2 <- c(vc2,v2)

        #derivability constraint
        v3 <- rep(0,nbT+1)
        v4 <- rep(0,nbT+1)
        v3[c(i,i+1)] <- c(deg * timeT[i]^(deg-1),-deg*timeT[i]^(deg-1))
        v4[c(i,i+1)] <- c((deg*timeT[i]^(deg-1))*sigma2[i],-(deg*timeT[i]^(deg-1))*sigma2[i+1])
        vd1 <- c(vd1,v3)
        vd2 <- c(vd2,v4)
        deg <- deg-1
      }
      m11[i,] <- c(vc1,rep(0,nbT*2))
      m12[i,] <- c(vd1,rep(0,nbT*2))
      m21[,i] <- vc2
      m22[,i] <- vd2
      i <- i+1
    }
    m1 <- rbind(m11,m12)
    m2 <- cbind(m21,m22)
  }

  # Models parameters
  m3 <- matrix(0,(nbT+1)*(d+1),(nbT+1)*(d+1))
  degre <- d
  m <- 1
  n <- 1

  while(m <= (nbT+1) & degre >=0){
    num <- m
    vp <- c()
    while(num <= (nbT+1)*(d+1)){
      vp <- c(vp,num)
      num <- num + (nbT+1)
    }
    deg2 <- d
    vect <- rep(0,(nbT+1)*(d+1))
    for(l in vp){
      vect[l] <- sum(X^(degre+deg2) * regimes_prob[,m])
      deg2 <- deg2-1
    }

    m3[n,] <- vect
    n <- n+1
    m <- m+1
    if(m > (nbT+1)){
      m <- 1
      degre <- degre-1
    }
  }

  if(constraint == 0){
    mFi <- m3
  }else{
    m4 <- cbind(m3,m2)
    mFi <- rbind(m1,m4)
  }

  return(mFi)

}
