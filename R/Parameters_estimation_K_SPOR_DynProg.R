
.Parameters_estimation_K_SPOR_DynProg <- function(Jacobian_M,vsolution,deg){

  param <- pseudoinverse(Jacobian_M, tol=sqrt(.Machine$double.eps)) %*% t(t(vsolution))

  nameMat = c()
  for(i in 1:(deg+1)){
    if((deg-i+1) == 0){
      nameMat = c(nameMat,"1")
    }else{
      nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
    }
  }

  mat_param <- matrix(param[1:(deg+1)],nrow=1,ncol=deg+1,dimnames = list(c(),nameMat))

  return(mat_param)
}
