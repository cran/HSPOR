
.Fixed_point_method_K_SPOR <- function(X,Y,deg,timeT,sigma2,obs_weights,regimes_prob,constraint,FP_nbIter=20){

  nbT <- length(timeT)
  for(p in 1:FP_nbIter){
    #vector parameters estimation
    M_mat <- .Jacobian_Matrix_K_SPOR(X,Y,deg,regimes_prob,timeT,sigma2,constraint)
    gammaV <- .Vector_solution(X,Y,deg,timeT,regimes_prob,constraint)
    vparam <- .Parameters_estimation(M_mat,gammaV)

    #Transformation to a matrix
    mat_param <- .Transform_vector_to_matrix(vparam,deg,timeT)

    #Estimation of sigma2
    sigma2 <- matrix(.Variance_estimation(X,Y,deg,timeT,mat_param,regimes_prob),nrow=1,ncol=nbT+1)
  }

  #Minus_Complete_log_likelihood
  MLL <- .Minus_Complete_log_likelihood(X,Y,deg,timeT,mat_param,sigma2,obs_weights)

  list(mat_param,sigma2,MLL)
}
