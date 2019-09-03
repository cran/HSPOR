

.Jacobian_Matrix_one_constraint_K_SPOR_DynProg <- function(X,Y,deg,sigma2,constraint,constraint_point){

  # Models parameters
  m3 <- matrix(0,deg+1,deg+1)
  d1 <- deg
  for(i in 1:(deg+1)){
    d2 <- deg
    for(j in 1:(deg+1)){
      m3[i,j] <- sum(X^(d1 + d2))
      d2 <- d2-1
    }
    d1 <- d1-1
  }

  #Constraints == 1 (continuity)
  if(constraint == 1){
    m1 <- matrix(0,1,(deg+1)+1)
    m2 <- matrix(0,(deg+1),1)

    for(i in 1:(deg+1)){
      m1[1,i] <- constraint_point[1]^(deg + 1 - i)
      m2[i,1] <- (constraint_point[1]^(deg + 1 - i))*sigma2
    }
  }

  #Constraints == 2 (continuity & derivability)
  if(constraint == 2){
    m1 <- matrix(0,2,(deg+1)+2)
    m2 <- matrix(0,(deg+1),2)

    for(i in 1:(deg+1)){
      m1[1,i] <- constraint_point[1]^(deg + 1 - i)
      m1[2,i] <- (deg + 1 - i)*(constraint_point[1]^(max(0,deg-i)))
      m2[i,1] <- (constraint_point[1]^(deg + 1 - i))*sigma2
      m2[i,2] <- (deg + 1 - i)*(constraint_point[1]^(max(0,deg-i)))*sigma2
    }
  }

  if(constraint == 0){
    matFi <- m3
  } else{
    m4 <- cbind(m3,m2)
    matFi <- rbind(m1,m4)
  }

  return(matFi)

}
