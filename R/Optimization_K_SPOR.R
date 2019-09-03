
.Optimization_K_SPOR <- function(X,Y,deg,timeT,EM,constraint){

  nbT <- length(timeT)

  #Parameters initialisation :
  mat_param_i <- .Init_param(X,Y,deg,timeT)[[1]]
  sigma2_i <- .Init_param(X,Y,deg,timeT)[[2]]

  #Optimisation
  if(constraint == 0){
    ML_estim <- .Maximum_likelihood_K_SPOR(X,Y,deg,timeT,mat_param_i,sigma2_i,constraint)
    param <- ML_estim[[1]]
    sigma2 <- ML_estim[[2]]
    MLL <- ML_estim[[3]]
  }else if((EM == FALSE) & (constraint >0)){
    obs_weights <- .Observations_weigths(X,timeT)
    regimes_prob <- .Regime_probability(X,Y,deg,timeT,mat_param_i,sigma2_i,obs_weights)
    FP_estim <- .Fixed_point_method_K_SPOR(X,Y,deg,timeT,sigma2_i,obs_weights,regimes_prob,constraint)
    param <- FP_estim[[1]]
    sigma2 <- FP_estim[[2]]
    MLL <- FP_estim[[3]]
  }else{
    EM_estim <- .EM_algorithm_K_SPOR(X,Y,deg,timeT,mat_param_i,sigma2_i,constraint)
    param <- EM_estim[[1]]
    sigma2 <- EM_estim[[2]]
    MLL <- EM_estim[[3]]
  }

  list(param,sigma2,MLL)
}
