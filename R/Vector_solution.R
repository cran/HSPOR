
.Vector_solution <- function(X,Y,deg,timeT,regimes_prob,constraint){

  nbT <- length(timeT)

  mat_X = cbind(sapply(deg:1,function(y) X^y),matrix(rep(1,length(X)),nrow = length(X),ncol=1))

  solV1 <-  matrix(t(regimes_prob)%*%(Y*mat_X),nrow <- 1,ncol <- (nbT+1)*(deg+1))

  if(constraint == 0){
    solV <- solV1
  }else if(constraint == 1){
    solV2 <- matrix(0,nrow=1,ncol = nbT)
    solV <- cbind(solV2,solV1)
  }else{
    solV3 <- matrix(0,nrow=1,ncol = nbT*2)
    solV <- cbind(solV3,solV1)
  }

  return(solV)
}
