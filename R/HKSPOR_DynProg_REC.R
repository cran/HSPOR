
.HKSPOR_DynProg_REC <- function(datX,datY,deg,EndPoint,K,constraint,grid,list_opt_K,smoothing_model,smoothing,TimeTX = c(),TimeTY=c(),parameters = data.frame(),Variances=c(),max_K = K){
  #Value of the K max
  max_K = max(max_K,K)

  #If results are already save, we keep these optimal results

  if(!is.null(list_opt_K[[K]][[which(grid == EndPoint[1])]]$MinusLogLikelihood)){
    best_estim <- c()
    sigma2 <- c()
    param <- c()
    MLL <- list_opt_K[[K]][[which(grid == EndPoint[1])]]$MinusLogLikelihood
    parameters <-  list_opt_K[[K]][[which(grid == EndPoint[1])]]$parameters
    Variances <- list_opt_K[[K]][[which(grid == EndPoint[1])]]$Variances
    TimeTX <- list_opt_K[[K]][[which(grid == EndPoint[1])]]$TimeTX
    TimeTY <- list_opt_K[[K]][[which(grid == EndPoint[1])]]$TimeTY

    #If we are studying only one regime or if we are in the last regime
  }else if(K == 1){

    if(smoothing==F){
      res = .Fixed_point_method_one_constraint_K_SPOR_DynProg(datX,datY,deg,1,constraint,EndPoint,'end')
    }else{
      res <- .Fixed_point_smoothing_method_one_constraint_K_SPOR_DynProg(datX,datY,deg,1,constraint,EndPoint,min(5,length(datX[which(datX<EndPoint[1])[length(which(datX<EndPoint[1]))]:length(datX)])-1),'end')
    }

    MLL <- res[[3]]
    best_estim <- c()
    sigma2 <- c()
    param <- c()
    parameters <- res[[1]]
    Variances <- res[[2]]
    TimeTX <- c()
    TimeTY <- c()

  }else{
    #x grid for the times of transition
    x_can = unique((datX[(((deg+1)*2*K)-(deg+1)):(length(which(datX < EndPoint[1]))-(deg+1))] + datX[(((deg+1)*2*K)-(deg+2)):(length(which(datX < EndPoint[1]))-(deg+2))])/2)

    #Y grid and yp grid associated to the x grid of times of transition
    if(length(x_can) > 1){
      y_can <- predict(smoothing_model,newdata = data.frame(x_can))$Estimation[,1]
      yp_can <- predict(smoothing_model,newdata = data.frame(x_can))$First_deriv[,1]
    }else{
      y_can <- predict(smoothing_model,newdata = data.frame(x_can))$Estimation[1]
      yp_can <- predict(smoothing_model,newdata = data.frame(x_can))$First_deriv[1]
    }

    opt = +Inf

    for (u in 1:length(x_can)){
      #print(u)
      trans_x = x_can[u]
      #print(length(datX[datX < trans_x]))
      trans_y = y_can[u]
      trans_yp <- yp_can[u]

      #If we work on the right regime : one constraint applied only at the beginning of the regime
      if((max_K == K) && (EndPoint[1] == datX[length(datX)])){
        if(smoothing == F){
          right_results <- .Fixed_point_method_one_constraint_K_SPOR_DynProg(datX,datY,deg,1,constraint,c(trans_x,trans_y,trans_yp),'beginning')
        }else{
          right_results <- .Fixed_point_smoothing_method_one_constraint_K_SPOR_DynProg(datX,datY,deg,1,constraint,c(trans_x,trans_y,trans_yp),min(5,length(datX[1:which(datX > trans_x)[1]])-1,length(datX[which(datX > trans_x)[1]:length(datX)])),'beginning')
        }

      }else{
        #For all the other regimes except the last one : two constraints applied at the beginning and at the end of the regime
        if(smoothing == F){
          right_results = .Fixed_point_method_two_constraints_K_SPOR_DynProg(datX,datY,deg,1,constraint,c(trans_x,trans_y,trans_yp),EndPoint)
        }else{
          right_results = .Fixed_point_smoothing_method_two_constraints_K_SPOR_DynProg(datX,datY,deg,1,constraint,c(trans_x,trans_y,trans_yp),EndPoint,
                                                                                       min(5,length(datX[1:which(datX>trans_x)[1]])-1,length(datX[which(datX<EndPoint[1])[length(which(datX<EndPoint[1]))]:length(datX)])-1
                                                                                           ,floor(length(datX[which(datX>trans_x & datX < EndPoint[1])])/2)))
        }
      }

      Right_MLL = right_results[[3]]

      # Dynamic programming
      left_results <- .HKSPOR_DynProg_REC(datX,datY,deg,c(trans_x,trans_y,trans_yp),K-1,constraint,grid,list_opt_K,smoothing_model,smoothing,TimeTX,TimeTY,parameters,Variances,max_K)
      MLL_gauche = left_results$MinusLogLikelihood

      MLL = Right_MLL + MLL_gauche
      #print(MLL)

      if(MLL < opt){
        best_estim = c(trans_x,trans_y)
        TimeTX <- left_results$TimeTX
        TimeTY <- left_results$TimeTY
        parameters <- left_results$parameters
        Variances <- left_results$Variances
        param <- right_results[[1]]
        sigma2 <- right_results[[2]]
        opt = MLL
      }
    }
    MLL = opt
  }


  #We save the results
  list_opt_K[[K]][[which(grid == EndPoint[1])]] <-  list("MinusLogLikelihood" = MLL, "parameters" = rbind(parameters,param), "Variances" = c(Variances,sigma2), "TimeTX" = c(TimeTX,best_estim[1]), "TimeTY" = c(TimeTY,best_estim[2]))

  #output
  list("MinusLogLikelihood" = MLL, "parameters" = rbind(parameters,param), "Variances" = c(Variances,sigma2), "TimeTX" = c(TimeTX,best_estim[1]), "TimeTY" = c(TimeTY,best_estim[2]),"list_opt_K" = list_opt_K)
}
