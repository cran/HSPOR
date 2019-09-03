

.Fixed_point_smoothing_method_two_constraints_K_SPOR_DynProg <- function(datX,datY,deg,sigma2,constraint,begin_point,end_point,nb_pt_pond,FP_nbIter=20){

  #Weight of the observations
  datX = unique(datX)
  if(nb_pt_pond != 0){
    X <- c(datX[(which(datX == datX[(datX < end_point[1]) & datX>=begin_point[1]][1])-nb_pt_pond):(which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][1])-1)]
           ,datX[datX<end_point[1] & datX>=begin_point[1]],
           datX[(which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][length(datX[datX<end_point[1] & datX>=begin_point[1]])])+1):
                  (which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][length(datX[datX<end_point[1] & datX>=begin_point[1]])])+nb_pt_pond)])

    Y <- c(datY[(which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][1])-nb_pt_pond):(which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][1])-1)]
           ,datY[datX<end_point[1] & datX>=begin_point[1]],
           datY[(which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][length(datX[datX<end_point[1] & datX>=begin_point[1]])])+1):
                  (which(datX == datX[datX<end_point[1] & datX>=begin_point[1]][length(datX[datX<end_point[1] & datX>=begin_point[1]])])+nb_pt_pond)])

    obs_weights1 <- 1-(1-((X[1:(2*nb_pt_pond)] - X[1])/(X[2*nb_pt_pond] - X[1]))^2)
    obs_weights2 <- 1-((X[max(1,(length(X)-(2*nb_pt_pond-1))):length(X)] - X[max(1,(length(X)-(2*nb_pt_pond-1)))])/(X[length(X)] - X[max(1,(length(X)-(2*nb_pt_pond-1)))]))^2
    obs_weights = c(obs_weights1,rep(1,(length(X)-length(c(obs_weights1,obs_weights2)))),obs_weights2)


  }else{
    X <- datX[datX < end_point[1] & datX>=begin_point[1]]
    Y <- datY[datX < end_point[1] & datX>=begin_point[1]]
    obs_weights <- rep(1,length(X))
  }

  for(p in 1:FP_nbIter){
    #vector parameters estimation
    M_mat <- .Jacobian_Matrix_two_constraints_K_SPOR_DynProg(X,Y,deg,sigma2,constraint,begin_point,end_point)
    gammaV <- .Vector_solution_two_constraints_K_SPOR_DynProg(X,Y,deg,constraint,begin_point,end_point)
    mat_param <- .Parameters_estimation_K_SPOR_DynProg(M_mat,gammaV,deg)

    sigma2 <- .Variance_estimation_K_SPOR_DynProg(X,Y,deg,mat_param)
  }

  #Complete_log_likelihood
  s <- 0
  for(i in 1:(deg+1)){
    s <- s + mat_param[1,i] * X^(deg+1-i)
  }

  wp <-  obs_weights*((1/sqrt(2*pi*sigma2)) * exp( - ((Y - s)^2)/(2*sigma2)))

  penalised_MLL <- -sum(log(wp[-c(1,length(wp))]))

  list(mat_param,sigma2,penalised_MLL)
}
