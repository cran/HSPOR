
.Jacobian_Matrix_many_constraints_2_SPOR_DynProg <- function(X,Y,deg,regimes_prob,timeT,sigma2,constraint,constraint_point,algo_direction){

  nbT <- length(timeT)
  d = deg

  if(constraint == 0){
    t_mat <- (nbT+1)*(d+1)
  }else if(constraint == 1){
    if(nrow(constraint_point)==1){
      t_mat <- (nbT+1)*(d+1) + (nbT+1)  #(nbT+1) parametres contraintes
    }else{
      t_mat <- (nbT+1)*(d+1) + (nbT+2)  #(nbT+2) parametres contraintes
    }
  }else{
    if(nrow(constraint_point)==1){
      t_mat <- (nbT+1)*(d+1) + 2*(nbT+1)  #2*(nbT+1) parametres contraintes
    }else{
      t_mat <- (nbT+1)*(d+1) + 2*(nbT+2)  #2*(nbT+2) parametres contraintes
    }
  }

  #Constraints
  if(constraint == 1){ #if continuity constraint

    if(nrow(constraint_point) == 1){
      m1 <- matrix(0,nbT+1,t_mat)
      m2 <- matrix(0,(nbT+1)*(d+1),nbT+1)

      v1 <- c()
      v2 <- c()
      v3 <- c()
      v4 <- c()

      for(i in 1 : (d+1)){
        v1 <- c(v1,c(timeT^(d + 1 - i),-timeT^(d + 1 - i)))
        v2 <- c(v2,c((timeT^(d + 1 - i))*sigma2[1],-(timeT^(d + 1 - i))*sigma2[2]))

        if(algo_direction == 'left'){
          v3 <- c(v3,c(0,constraint_point[1,1]^(d + 1 - i)))
          v4 <- c(v4,c(0,(constraint_point[1,1]^(d + 1 - i))*sigma2[2]))
        }else{
          v3 <- c(v3,c(constraint_point[1,1]^(d + 1 - i),0))
          v4 <- c(v4,c((constraint_point[1,1]^(d + 1 - i))*sigma2[1],0))
        }
      }

      m1[1,] <- c(v1,rep(0,nbT+1))
      m1[2,] <- c(v3,rep(0,nbT+1))
      m2[,1] <- v2
      m2[,2] <- v4

    }else{
      m1 <- matrix(0,nbT+2,t_mat)
      m2 <- matrix(0,(nbT+1)*(d+1),nbT+2)

      v1 <- c()
      v2 <- c()
      v3 <- c()
      v4 <- c()
      v5 <- c()
      v6 <- c()

      for(i in 1:(d+1)){

        v1 <- c(v1,c(timeT^(d+1-i),-timeT^(d+1-i)))
        v2 <- c(v2,c((timeT^(d+1-i))*sigma2[1],-(timeT^(d+1-i))*sigma2[2]))

        v3 <- c(v3,c(constraint_point[1,1]^(d+1-i),0))
        v4 <- c(v4,c((constraint_point[1,1]^(d+1-i))*sigma2[1],0))

        v5 <- c(v5,c(0,constraint_point[2,1]^(d+1-i)))
        v6 <- c(v6,c(0,(constraint_point[2,1]^(d+1-i))*sigma2[2]))
      }

      m1[1,] <- c(v1,rep(0,nbT+2))
      m1[2,] <- c(v3,rep(0,nbT+2))
      m1[3,] <- c(v5,rep(0,nbT+2))
      m2[,1] <- v2
      m2[,2] <- v4
      m2[,3] <- v6
    }
  }


  if(constraint == 2){ #if continuity constraint

    if(nrow(constraint_point) == 1){
      m1 <- matrix(0,2*(nbT+1),t_mat)
      m2 <- matrix(0,(nbT+1)*(d+1),2*(nbT+1))

      v1c <- c()
      v2c <- c()
      v3c <- c()
      v4c <- c()

      v1d <- c()
      v2d <- c()
      v3d <- c()
      v4d <- c()

      for(i in 1:(d+1)){
        v1c <- c(v1c,c(timeT^(d+1-i),-timeT^(d+1-i)))
        v2c <- c(v2c,c((timeT^(d+1-i))*sigma2[1],-(timeT^(d+1-i))*sigma2[2]))
        v1d <- c(v1d,c((d+1-i)*timeT^(d-i),-(d+1-i)*timeT^(d-i)))
        v2d <- c(v2d,c(((d+1-i)*timeT^(d-i))*sigma2[1],-((d+1-i)*timeT^(d-i))*sigma2[2]))

        if(algo_direction == 'left'){
          v3c <- c(v3c,c(0,constraint_point[1,1]^(d+1-i)))
          v4c <- c(v4c,c(0,(constraint_point[1,1]^(d+1-i))*sigma2[2]))

          v3d <- c(v3d,c(0,(d+1-i)*constraint_point[1,1]^(d-i)))
          v4d <- c(v4d,c(0,((d+1-i)*constraint_point[1,1]^(d-i))*sigma2[2]))
        }else{
          v3c <- c(v3c,c(constraint_point[1,1]^(d+1-i),0))
          v4c <- c(v4c,c((constraint_point[1,1]^(d+1-i))*sigma2[1],0))

          v3d <- c(v3d,c((d+1-i)*constraint_point[1,1]^(d-i),0))
          v4d <- c(v4d,c(((d+1-i)*constraint_point[1,1]^(d-i))*sigma2[1],0))
        }
      }

      m1[1,] <- c(v1c,rep(0,2*(nbT+1)))
      m1[2,] <- c(v3c,rep(0,2*(nbT+1)))
      m1[3,] <- c(v1d,rep(0,2*(nbT+1)))
      m1[4,] <- c(v3d,rep(0,2*(nbT+1)))
      m2[,1] <- v2c
      m2[,2] <- v4c
      m2[,4] <- v4d

    }else{
      m1 <- matrix(0,2*(nbT+2),t_mat)
      m2 <- matrix(0,(nbT+1)*(d+1),2*(nbT+2))

      v1c <- c()
      v2c <- c()
      v3c <- c()
      v4c <- c()
      v5c <- c()
      v6c <- c()

      v1d <- c()
      v2d <- c()
      v3d <- c()
      v4d <- c()
      v5d <- c()
      v6d <- c()

      for(i in 1:(d+1)){
        #continuity between the two regimes
        v1c <- c(v1c,c(timeT^(d+1-i),-timeT^(d+1-i)))
        v2c <- c(v2c,c((timeT^(d+1-i))*sigma2[1],-(timeT^(d+1-i))*sigma2[2]))
        #continuite ? gauche
        v3c <- c(v3c,c(constraint_point[1,1]^(d+1-i),0))
        v4c <- c(v4c,c((constraint_point[1,1]^(d+1-i))*sigma2[1],0))
        #continuite ? droite
        v5c <- c(v5c,c(0,constraint_point[2,1]^(d+1-i)))
        v6c <- c(v6c,c(0,(constraint_point[2,1]^(d+1-i))*sigma2[2]))

        #differentiability constraints
        v1d <- c(v1d,c((d+1-i)*timeT^(d-i),-(d+1-i)*timeT^(d-i)))
        v2d <- c(v2d,c(((d+1-i)*timeT^(d-i))*sigma2[1],-((d+1-i)*timeT^(d-i))*sigma2[2]))
        #differentiability constrainton the left
        v3d <- c(v3d,c((d+1-i)*constraint_point[1,1]^(d-i),0))
        v4d <- c(v4d,c(((d+1-i)*constraint_point[1,1]^(d-i))*sigma2[1],0))
        #differentiability constraint on the right
        v5d <- c(v5d,c(0,(d+1-i)*constraint_point[2,1]^(d-i)))
        v6d <- c(v6d,c(0,((d+1-i)*constraint_point[2,1]^(d-i))*sigma2[2]))
      }

      m1[1,] <- c(v1c,rep(0,2*(nbT+2)))
      m1[2,] <- c(v3c,rep(0,2*(nbT+2)))
      m1[3,] <- c(v5c,rep(0,2*(nbT+2)))
      m1[4,] <- c(v1d,rep(0,2*(nbT+2)))
      m1[5,] <- c(v3d,rep(0,2*(nbT+2)))
      m1[6,] <- c(v5d,rep(0,2*(nbT+2)))
      m2[,1] <- v2c
      m2[,2] <- v4c
      m2[,3] <- v6c
      m2[,4] <- v2d
      m2[,5] <- v4d
      m2[,6] <- v6d
    }
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
    JacobianMat <- m3
  }else{
    m4 <- cbind(m3,m2)
    JacobianMat <- rbind(m1,m4)
  }

  return(JacobianMat)

}
