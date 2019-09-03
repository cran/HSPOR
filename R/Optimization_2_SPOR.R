
.Optimization_2_SPOR <- function(X,Y,deg,timeT,mat_param,sigma2,constraint,EM){

  if(constraint == 1 | constraint == 2){
    if(EM==TRUE){
      param_est <- .EM_algorithm_2_SPOR(X,Y,deg,timeT,mat_param,sigma2,constraint)
    }else{
      obs_weights_i <- .Observations_weigths(X,timeT)
      regimes_prob_i <- obs_weights_i
      param_est <- .Fixed_point_method_2_SPOR(X,Y,deg,timeT,sigma2,obs_weights_i,regimes_prob_i,constraint)
    }
  }else{
    obs_weights_i <- .Observations_weigths(X,timeT)
    MLL <- .Minus_Complete_log_likelihood(X,Y,deg,timeT,mat_param,sigma2,obs_weights_i)
    param_est <- list(mat_param,sigma2,MLL)
  }
  return(param_est)
}
