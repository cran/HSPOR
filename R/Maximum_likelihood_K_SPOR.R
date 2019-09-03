
.Maximum_likelihood_K_SPOR <- function(X,Y,deg,timeT,mat_param,sigma2,constraint){

  nbT <- length(timeT)

  obs_weights <- .Observations_weigths(X,timeT)
  regimes_prob <- .Regime_probability(X,Y,deg,timeT,mat_param,sigma2,obs_weights)

  M_mat <- .Jacobian_Matrix_K_SPOR(X,Y,deg,regimes_prob,timeT,sigma2,constraint)
  gammaV <- .Vector_solution(X,Y,deg,timeT,regimes_prob,constraint)
  vparam <- .Parameters_estimation(M_mat,gammaV)

  #Transformation to a matrix
  mat_param <- .Transform_vector_to_matrix(vparam,deg,timeT)

  #Estimation of sigma2
  sigma2 <- matrix(.Variance_estimation(X,Y,deg,timeT,mat_param,regimes_prob),nrow=1,ncol=nbT+1)

  #Minus complete log likelihood
  MLL <- .Minus_Complete_log_likelihood(X,Y,deg,timeT,mat_param,sigma2,obs_weights)

  list(mat_param,sigma2,MLL)
}
