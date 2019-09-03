
.Trichotomic_algorithm_many_constraints_2_SPOR_DynProg <- function(X,Y,deg,constraint,EM,constraint_point,algo_direction,stop_tricho=0.001){

  if(length(X) <= ((deg+2)*4)){
    intvalle = data.frame(t(c(X[deg+2],X[length(X)-(deg+2)])))
  }else{
    intvalle = rbind(c(X[deg+2],X[length(X)-(deg+2)]),c(X[floor(length(X)/2)],X[length(X)-(deg+2)]),
                     c(X[deg+2], X[floor(length(X)/2)]))
  }

  #MLL_Ti <- c()
  Topt <- c()
  mat_param <- list()
  sigma2 <- list()
  MLL <- c()


  for(m in 1 : length(intvalle[,1])){
    #Jumps initialization
    h = (intvalle[m,2]-intvalle[m,1])/3
    T1 = intvalle[m,1] + h
    T2 = intvalle[m,2] - h

    mat_parami_T1 <- .Init_param(X,Y,deg,T1)[[1]]
    sigma2i_T1 <- .Init_param(X,Y,deg,T1)[[2]]

    mat_parami_T2 <- .Init_param(X,Y,deg,T2)[[1]]
    sigma2i_T2 <- .Init_param(X,Y,deg,T2)[[2]]


    while(intvalle[m,2]-intvalle[m,1] > stop_tricho){
      opt_T1 <- .Optimization_many_constraints_2_SPOR_DynProg(X,Y,deg,T1,mat_parami_T1,sigma2i_T1,EM,constraint,constraint_point,algo_direction)
      opt_T2 <- .Optimization_many_constraints_2_SPOR_DynProg(X,Y,deg,T2,mat_parami_T2,sigma2i_T2,EM,constraint,constraint_point,algo_direction)

      #MLL_Ti <- c(MLL_Ti,opt_T1[[3]],opt_T2[[3]])

      if(opt_T1[[3]] <= opt_T2[[3]]){
        intvalle[m,2] <- T2
        Topt[m] = T1
        mat_param[[m]] <- opt_T1[[1]]
        sigma2[[m]] <- opt_T1[[2]]
        MLL[m] <- opt_T1[[3]]
        #MLL_Ti <- c(MLL_Ti,MLL[m])
      }
      if(opt_T1[[3]] >= opt_T2[[3]]){
        intvalle[m,1]= T1
        Topt[m] = T2
        mat_param[[m]] = opt_T2[[1]]
        sigma2[[m]] = opt_T2[[2]]
        MLL[m] <- opt_T2[[3]]
      }

      h = (intvalle[m,2] - intvalle[m,1])/3
      T1 = intvalle[m,1] + h
      T2 = intvalle[m,2] - h

      mat_parami_T1 <- .Init_param(X,Y,deg,T1)[[1]]
      sigma2i_T1 <- .Init_param(X,Y,deg,T1)[[2]]

      mat_parami_T2 <- .Init_param(X,Y,deg,T2)[[1]]
      sigma2i_T2 <- .Init_param(X,Y,deg,T2)[[2]]

    }
  }

  #dispersion_MLL <- sqrt(var(MLL_Ti))

  list(Topt[which.min(MLL)],mat_param[[which.min(MLL)]],sigma2[[which.min(MLL)]],MLL[ which.min(MLL)])
}
