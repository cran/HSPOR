
.Vector_solution_one_constraint_2_SPOR_DynProg <- function(X,Y,deg,timeT,regimes_prob,constraint){

  nbT <- length(timeT)
  #ıCalcul du vecteur sol
  mat_X = cbind(sapply(deg:1,function(y) X^y),matrix(rep(1,length(X)),nrow = length(X),ncol=1))

  solV1 <-  matrix(t(regimes_prob)%*%(Y*mat_X),nrow <- 1,ncol <- (nbT+1)*(deg+1))

  if(constraint == 0){
    solV <- solV1
  }else if(constraint == 1){
    solV2 <- matrix(0,nrow=1,ncol = nbT)
    solV <- cbind(solV2,solV1)
  }else{
    #continuity
    solV2 <- matrix(c(0,0),nrow=1,ncol = nbT+1)
    solV <- cbind(solV2,solV1)
  }

  return(solV)
}
