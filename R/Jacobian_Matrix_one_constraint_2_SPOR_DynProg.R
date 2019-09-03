

.Jacobian_Matrix_one_constraint_2_SPOR_DynProg <- function(X,Y,deg,regimes_prob,timeT,sigma2,constraint){

  nbT <- length(timeT)
  d = deg

  if(constraint == 0){
    t_mat <- (nbT+1)*(d+1)
  }else if(constraint == 1){
    t_mat <- (nbT+1)*(d+1) + nbT  #(nbT) param?tres contraintes
  }else{
    t_mat <- (nbT+1)*(d+1) + (2*nbT) #(nbT)*2 param?tres de contraintes
  }

  #Constraints
  if(constraint == 1){ #if continuity constraint
    m1 <- matrix(0,nbT,t_mat)
    m2 <- matrix(0,(nbT+1)*(d+1),nbT)

    v1 <- c()
    v2 <- c()
    for(i in 1:(d+1)){
      v1 <- c(v1,c(timeT^(d+1-i),-timeT^(d+1-i)))
      v2 <- c(v2,c((timeT^(d+1-i))*sigma2[1],-(timeT^(d+1-i))*sigma2[2]))
    }
    m1[1,] <- c(v1,0)
    m2[,1] <- v2
  }

  if(constraint == 2){ #if continuity constraint
    m11 <- matrix(0,nbT,t_mat)
    m12 <- matrix(0,nbT,t_mat)
    m21 <- matrix(0,(nbT+1)*(d+1),nbT)
    m22 <- matrix(0,(nbT+1)*(d+1),nbT)

    vc1 <- c()
    vc2 <- c()
    vd1 <- c()
    vd2 <- c()

    for(i in 1:(d+1)){
      vc1 <- c(vc1,c(timeT^(d+1-i),-timeT^(d+1-i)))
      vc2 <- c(vc2,c((timeT^(d+1-i))*sigma2[1],-(timeT^(d+1-i))*sigma2[2]))

      vd1 <- c(vd1,c((d+1-i) * timeT^(d-i),-(d+1-i)*timeT^(d-i)))
      vd2 <- c(vd2,c(((d+1-i)*timeT^(d-i))*sigma2[1],-((d+1-i)*timeT^(d-i))*sigma2[2]))
    }

    m11[1,] <- c(vc1,rep(0,2))
    m21[,1] <- vc2
    m12[1,] <- c(vd1,rep(0,2))
    m22[,1] <- vd2

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
