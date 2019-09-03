
.EM_algorithm_2_SPOR <- function(X,Y,deg,timeT,mat_param,sigma2,constraint,EM_nbIter=20){

  for(j in 1:EM_nbIter){

    obs_weights <- .Smoothing_observations_weigths_2_SPOR(X,deg,timeT)
    regimes_prob <- .Regime_probability(X,Y,deg,timeT,mat_param,sigma2,obs_weights)

    FP_estimation <- .Fixed_point_method_2_SPOR(X,Y,deg,timeT,sigma2,obs_weights,regimes_prob,constraint)

    mat_param <- FP_estimation[[1]]
    sigma2 <- FP_estimation[[2]]

  }

  MLL <- .Minus_Complete_log_likelihood(X,Y,deg,timeT,mat_param,sigma2,obs_weights)

  list(mat_param,sigma2,MLL)
}
