
.Vector_solution_many_constraints_2_SPOR_DynProg <- function(X,Y,deg,TimeT,regimes_prob,constraint,constraint_point){

  nbT <- length(TimeT)
  #solution for the parameters
  mat_X = cbind(sapply(deg:1,function(y) X^y),matrix(rep(1,length(X)),nrow = length(X),ncol=1))

  solV1 <-  matrix(t(regimes_prob)%*%(Y*mat_X),nrow <- 1,ncol <- (nbT+1)*(deg+1))

  #solution for the constraints
  if(constraint == 0){
    solV <- solV1
  }else if(constraint == 1){
    if( nrow(constraint_point) == 1 ){
      solV2 <- matrix(c(0,constraint_point[1,2]),nrow=1,ncol=nbT+1)
    }else{
      solV2 <- matrix(c(0,constraint_point[1,2],constraint_point[2,2]),nrow=1,ncol=nbT+2)
    }
    solV <- cbind(solV2,solV1)
  }else{
    if( nrow(constraint_point) == 1 ){
      solV2 <- matrix(c(0,constraint_point[1,2]),nrow=1,ncol=nbT+1)
      solV3 <- matrix(c(0,constraint_point[1,3]),nrow=1,ncol=nbT+1)
    }else{
      solV2 <- matrix(c(0,constraint_point[1,2],constraint_point[2,2]),nrow=1,ncol=nbT+2)
      solV3 <- matrix(c(0,constraint_point[1,3],constraint_point[2,3]),nrow=1,ncol=nbT+3)
    }
    solV <- cbind(solV2,solV3,solV1)
  }

  return(solV)
}
