
.Fixed_point_smoothing_method_one_constraint_K_SPOR_DynProg <- function(datX,datY,deg,sigma2,constraint,constraint_point,nb_pt_pond,indic = 'beginning',FP_nbIter=20){
  if(indic == 'beginning'){
    if(nb_pt_pond != 0){
      #lissage autour du saut : on prend 5 observations supplémentaires avant l'instant de transition
      X <- c(datX[(which(datX == datX[datX>=constraint_point[1]][1])-nb_pt_pond):(which(datX == datX[datX>=constraint_point[1]][1])-1)],datX[datX>=constraint_point[1]])
      Y <- c(datY[(which(datX == datX[datX>=constraint_point[1]][1])-nb_pt_pond):(which(datX == datX[datX>=constraint_point[1]][1])-1)],datY[datX>=constraint_point[1]])

      obs_weights1 <- 1-(1-((X[1:(2*nb_pt_pond)] - X[1])/(X[2*nb_pt_pond] - X[1]))^2)
      obs_weights = c(obs_weights1,rep(1,length(X)-length(obs_weights1)))

    }else{
      X <- datX[datX>=constraint_point[1]]
      Y <- datY[datX>=constraint_point[1]]
      obs_weights <- rep(1,length(X))
    }

  }else {

    #ponderation de la vraisemblance
    if(nb_pt_pond != 0){
      #Si point de debut on prend 5 observations supplementaires après l'instant de transition
      X <- c(datX[datX<constraint_point[1]],datX[(which(datX == datX[datX<constraint_point[1]][length(datX[datX<constraint_point[1]])])+1):
                                              (which(datX == datX[datX<constraint_point[1]][length(datX[datX<constraint_point[1]])])+nb_pt_pond)])
      Y <- c(datY[datX<constraint_point[1]],datY[(which(datX == datX[datX<constraint_point[1]][length(datX[datX<constraint_point[1]])])+1):
                                              (which(datX == datX[datX<constraint_point[1]][length(datX[datX<constraint_point[1]])])+nb_pt_pond)])

      obs_weights1 <- 1-((X[max(1,(length(X)-(2*nb_pt_pond-1))):length(X)] - X[max(1,(length(X)-(2*nb_pt_pond-1)))])/(X[length(X)] - X[max(1,(length(X)-(2*nb_pt_pond-1)))]))^2
      obs_weights = c(rep(1,length(X)-length(obs_weights1)),obs_weights1)
    }else{
      X <- datX[datX<constraint_point[1]]
      Y <- datY[datX<constraint_point[1]]
      obs_weights <- rep(1,length(X))
    }
  }

  for(p in 1:FP_nbIter){
    # parameters estimation
    M_mat <- .Jacobian_Matrix_one_constraint_K_SPOR_DynProg(X,Y,deg,sigma2,constraint,constraint_point)
    gammaV <- .Vector_solution_one_constraint_K_SPOR_DynProg(X,Y,deg,constraint,constraint_point)
    mat_param <- .Parameters_estimation_K_SPOR_DynProg(M_mat,gammaV,deg)

    sigma2 <- .Variance_estimation_K_SPOR_DynProg(X,Y,deg,mat_param)
  }

  #penalized Minus_Complete_log_likelihood
  s <- 0
  for(i in 1:(deg+1)){
    s <- s + mat_param[1,i] * X^(deg+1-i)
  }

  wp <-  obs_weights*((1/sqrt(2*pi*sigma2)) * exp( - ((Y - s)^2)/(2*sigma2)))

  if (indic == 'beginning'){
    penalised_MLL <- -sum(log(wp[-1]))
  }else{
    penalised_MLL <- -sum(log(wp[-length(wp)]))
  }

  list(mat_param,sigma2,penalised_MLL)
}
