
.Parameters_estimation <- function(Jacobian_M,vsolution){

  # Parameters estimation
  param <- pseudoinverse(Jacobian_M, tol=sqrt(.Machine$double.eps)) %*% as.vector(vsolution)

  return(param)
}
