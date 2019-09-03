
.Trichotomic_algorithm_2_SPOR <- function(X,Y,deg,constraint,EM,stop_tricho=0.001){

  intvalle = c(X[deg+2],X[length(X)-(deg+2)])

  #Jumps initialisations
  h = (intvalle[2]-intvalle[1])/3
  T1 = intvalle[1] + h
  T2 = intvalle[2] - h

  mat_parami_T1 <- .Init_param(X,Y,deg,T1)[[1]]
  sigma2i_T1 <- .Init_param(X,Y,deg,T1)[[2]]

  mat_parami_T2 <- .Init_param(X,Y,deg,T2)[[1]]
  sigma2i_T2 <- .Init_param(X,Y,deg,T2)[[2]]

  while(intvalle[2]-intvalle[1] > stop_tricho){
    opt_T1 <- .Optimization_2_SPOR(X,Y,deg,T1,mat_parami_T1,sigma2i_T1,constraint,EM)
    opt_T2 <- .Optimization_2_SPOR(X,Y,deg,T2,mat_parami_T2,sigma2i_T2,constraint,EM)

    if(opt_T1[[3]] <= opt_T2[[3]]){
      intvalle[2] <- T2
      Topt = T1
      mat_param <- opt_T1[[1]]
      sigma2 <- opt_T1[[2]]
      MLL <- opt_T1[[3]]
    }
    if(opt_T1[[3]] >= opt_T2[[3]]){
      intvalle[1]= T1
      Topt = T2
      mat_param = opt_T2[[1]]
      sigma2 = opt_T2[[2]]
      MLL <- opt_T2[[3]]
    }

    h = (intvalle[2]-intvalle[1])/3
    T1 = intvalle[1] + h
    T2 = intvalle[2] - h

    mat_parami_T1 <- .Init_param(X,Y,deg,T1)[[1]]
    sigma2i_T1 <- .Init_param(X,Y,deg,T1)[[2]]

    mat_parami_T2 <- .Init_param(X,Y,deg,T2)[[1]]
    sigma2i_T2 <- .Init_param(X,Y,deg,T2)[[2]]

  }

  list(Topt,mat_param,sigma2,MLL)
}
