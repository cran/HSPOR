

.Fixed_point_method_two_constraints_K_SPOR_DynProg <- function(datX,datY,deg,sigma2,constraint,begin_point,end_point,FP_nbIter=20){

  X = datX[datX<=end_point[1] & datX>begin_point[1]]
  Y = datY[datX<=end_point[1] & datX>begin_point[1]]

  for(p in 1:FP_nbIter){
    #vector parameters estimation
    M_mat <- .Jacobian_Matrix_two_constraints_K_SPOR_DynProg(X,Y,deg,sigma2,constraint,begin_point,end_point)
    gammaV <- .Vector_solution_two_constraints_K_SPOR_DynProg(X,Y,deg,constraint,begin_point,end_point)
    mat_param <- .Parameters_estimation_K_SPOR_DynProg(M_mat,gammaV,deg)

    sigma2 <- .Variance_estimation_K_SPOR_DynProg(X,Y,deg,mat_param)
  }

  #Minus_Complete_log_likelihood
  s <- 0
  for(i in 1:(deg+1)){
    s <- s + mat_param[1,i] * X^(deg+1-i)
  }

  pZ = rep(1,length(X))
  wp <-  pZ*((1/sqrt(2*pi*sigma2)) * exp( - ((Y - s)^2)/(2*sigma2)))

  MLL <- -sum(log(wp))

  list(mat_param,sigma2,MLL)
}
