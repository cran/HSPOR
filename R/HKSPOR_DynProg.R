#' @title Inference method for any number K of regimes using dynamic programming
#'
#' @description HKSPOR_DynProg is an inference method implemented in the form of a Bellman algorithm that estimates, under the assumption of regularity,
#' the parameters of a polynomial regression model for a number K of regimes given by the user..
#'
#' @param X A numerical vector corresponding to the explanatory variable. X must be sorted in ascending order
#' if this is not the case, X will be sorted in the function and the corresponding permutation will be applied to Y. The
#' user will be notified by a warning message. In addition, if X contains NAs, they will be deleted from the data and the user will be notified by a warning message.
#' Finally, if X contains duplicate data, the excess data will be deleted and the user will be notified by a warning message.
#' @param Y A numerical vector corresponding to the variable to be explain. It should contain at least two regimes that
#' could be modelled by polynomials. In addition, if Y contains NAs they will be deleted from the data and the
#' user will be notified by a warning message. Finally, if X contains dupplicate data, the excess data will be deleted and
#' the value of the remaining Y will become the average of the Ys, calculated for this value of X.
#' @param deg The degree of the polynomials. The size of X and Y must be greater than K(deg+2) + K.
#' @param K The number of regimes. The size of X and Y must be greater than K(deg+2) + K.
#' @param constraint  Number that determines the regularity assumption that is applied for the parameters estimation.
#' By default, the variable is set to 1, i. e. the parameters estimation is done under continuity constraint.
#' If the variable is 0 or 2, the estimation of the parameters will be done without assumption of regularity
#' (constraint = 0) or under assumption of differentiability (constraint = 2). Warning, if the differentiability
#' assumption is not verified by the model, it is preferable not to use it to estimate the model parameters.
#' In addition, in this dynamic programming method, to ensure that the number of constraints is not greater
#' that the number of parameters to be estimated, the degree of the polynomials must be at least equal to
#' 3 to be able to use the differentiability assumption.
#' @param smoothing A Boolean. If TRUE (default), the method will estimate the parameters of a piecewise polynomial regression
#' model with latent variable by maximizing the log-likelihood weighted by the probability of being in the
#' latent variable regime. If FALSE, the method will estimate the parameters of the piecewise polynomial regression
#' model.
#' @param verbose A Boolean. If FALSE (default) the HKSPOR_Dynprog function will return only one dataframe
#' containing the parameter estimates obtained for a model at K regimes. If TRUE, the function will return
#' all the results obtained for a model with 1 regime up to K regimes.
#' @param plotG A Boolean. If TRUE (default) the estimation results obtained by the HKSPOR_DynProg function are plotted.
#'
#' @return One or more dataframes depend on the verbose value. If verbose = False, the output table will
#' contain the estimated parameters of the polynomial regression model at K regimes: jump times, polynomial
#' coefficients and variances of K regimes. If verbose = True then there will be K dataframes in output.
#' Each table will contain the results of the estimated parameters obtained for each value of k (k=1,...,k=K).
#' If plotG = TRUE, the data (X,Y) and the estimated model(s) will
#' be plotted.
#'
#' @import npregfast
#' @import graphics
#' @import corpcor
#' @import stats
#' @export
#'
#' @examples
#' #generated data with three regimes
#' set.seed(1)
#' xgrid1 = seq(0,10,length.out=6)
#' xgrid2 = seq(10.2,20,length.out=6)

#' ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,4)
#' ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,4)

#' datX = c(xgrid1,xgrid2)
#' datY = c(ygrid1,ygrid2)
#'
#' #Inference of a polynomial regression model with two regimes (K=2) on these data.
#' #The degree of the polynomials is fixed to 2 and the parameters are estimated
#' #under continuity constraint.
#' HKSPOR_DynProg(datX,datY,2,2)
#'
#' \donttest{
#' set.seed(2)
#' xgrid1 = seq(0,10,by=0.2)
#' xgrid2 = seq(10.2,20,by=0.2)
#' xgrid3 = seq(20.2,30,by=0.2)

#' ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,3)
#' ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,3)
#' ygrid3 = -10*xgrid3+300+rnorm(length(xgrid3),0,3)

#' datX = c(xgrid1,xgrid2,xgrid3)
#' datY = c(ygrid1,ygrid2,ygrid3)
#'
#' #Inference of a polynomial regression model with three (K=3) regimes on these data.
#' #The degree of the polynomials is fixed to 2 and the parameters are estimated
#' #under continuity constraint.
#' HKSPOR_DynProg(datX,datY,2,3)
#' #Executed time : 3.658121 mins (intel core i7 processor)
#' }


HKSPOR_DynProg <- function(X,Y,deg,K,constraint=1,smoothing=TRUE,verbose=FALSE,plotG=TRUE){

  if( (length(which(is.na(X))) != 0)){
    Y = Y[-which(is.na(X))]
    X = X[-which(is.na(X))]
    warning("X contains missing data, these data have been deleted")
  }
  if( (length(which(is.na(Y))) != 0)){
    X = X[-which(is.na(Y))]
    Y = Y[-which(is.na(Y))]
    warning("Y contains missing data, these data have been deleted")
  }
  if(length(which(X != X[order(X)])) > 0){
    Y = Y[order(X)]
    X = X[order(X)]
    warning("X will be sorted and the corresponding permutation will be applied to Y")
  }

  if(length(X) != length(unique(X))){
    duble_value = unique(X[which(duplicated(X)==TRUE)])
    for(i in 1:length(duble_value)){
      new_y = sum(Y[which(X == duble_value[i])])/length(which(X == duble_value[i]))
      Y[which(X == duble_value[i])[1]] = new_y
      Y = Y[-which(X == duble_value[i])[2:length(which(X == duble_value[i]))]]
      X = X[-which(X == duble_value[i])[2:length(which(X == duble_value[i]))]]
      warning("Duplicate data of X have been delated. Y becames the mean of these data")
    }
  }

  if(length(X) < K*(deg+2) + K){
    warning(paste("X must contain at least",K*(deg+2)+K,"observations to use this degree and this number of regimes"))
    deg = floor((length(X)-(3*K))/K)
    if(deg <= 0){
      deg = 1
      K = floor(length(X)/(deg+3))
      warning(paste("The degree was reduced to 1 and the number of regimes to",K))
    }
    warning(paste("The degree was reduced to", deg))
  }

  if(length(X) < 3*K + K ){
    stop(paste("X must contain at least", 3*K + K, "observations for a degree equal to 1 and a number", K, "of regime"))
  }

  list_res <- list()
  x_grid <- c(X,sapply(1:(length(X)-1), function(i) (X[i]+X[i+1])/2))
  x_grid <- unique(x_grid[order(x_grid)])

  smoothing_model <- frfast(Y ~ X, p = deg, seed = 130853)

  #Estimation of Y and Y' at the end point by a smoothing model
  endPoint = X[length(X)]
  y_endPoint <- predict(smoothing_model,newdata = data.frame(endPoint))$Estimation[1]
  yp_endPoint <- predict(smoothing_model,newdata = data.frame(endPoint))$First_deriv[1]

  endPoint <- c(endPoint,y_endPoint,yp_endPoint)

  #Double list in which we could save the results for one K and and one end point
  list_opt_K <- list()
  for(l in 1 : K){
    list_opt_K[[l]] <- list()
    for(j in 1:length(x_grid)){
      list_opt_K[[l]][[j]] <- list("MinusLogLikelihood","parameters","Variances","TimeTX","TimeTY")
    }
  }

  #For many end points and a value of K we run the algorithm
  for(i in 1:K){
    #print(paste("i",i,sep=" "))
    if(((deg+1)*2*i) < (length(which(X <= endPoint[1])))){
      x_grid_endPoint = unique((X[((deg+1)*2*i):(length(which(X <= endPoint[1])))] + X[(((deg+1)*2*i)-1):(length(which(X <= endPoint[1]))-1)])/2)
      if(length(x_grid_endPoint) > 1){
        y_grid_endPoint <- predict(smoothing_model,newdata = data.frame(x_grid_endPoint))$Estimation[,1]
        yp_grid_endPoint <- predict(smoothing_model,newdata = data.frame(x_grid_endPoint))$First_deriv[,1]
      }else{
        y_grid_endPoint <- predict(smoothing_model,newdata = data.frame(x_grid_endPoint))$Estimation[1]
        yp_grid_endPoint <- predict(smoothing_model,newdata = data.frame(x_grid_endPoint))$First_deriv[1]
      }
      for(l in 1:length(x_grid_endPoint)){
        #print(paste("l",l,sep = " "))
        result <- .HKSPOR_DynProg_REC(X,Y,deg,c(x_grid_endPoint[l],y_grid_endPoint[l],yp_grid_endPoint[l]),i,constraint,x_grid,list_opt_K,smoothing_model,smoothing)
        list_opt_K <- result$list_opt_K
      }
    }
    list_res_int <- .HKSPOR_DynProg_REC(X,Y,deg,endPoint,i,constraint,x_grid,list_opt_K,smoothing_model,smoothing)
    list_opt_K <- list_res_int$list_opt_K
    list_res[[i]] <- list("MinusLogLikelihood"=list_res_int$MinusLogLikelihood,"parameters"=list_res_int$parameters,"Variances"=list_res_int$Variances,"TimeTX"=list_res_int$TimeTX)
  }

  if(verbose == TRUE){

    nameMat = c()
    for(i in 1:(deg+1)){
      if((deg-i+1) == 0){
        nameMat = c(nameMat,"1")
      }else{
        nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
      }
    }

    sortie = list()
    for( l in 1 : length(list_res)){
      colnames(list_res[[l]]$parameters) = nameMat

      sortie[[l]] = cbind.data.frame(data.frame("TimeT" = c(X[1],list_res[[l]]$TimeTX )),
                                     list_res[[l]]$parameters,"sigma2" = list_res[[l]]$Variances,
                                     row.names=as.factor(c(1:(length(list_res[[l]]$TimeTX )+1))))

    }

    if(plotG == TRUE){

      for(j in 1:length(sortie)){
        plot(X,Y,pch=20,cex.main = 3, cex.lab = 2,cex.axis = 2)
        if(length(sortie[[j]]$TimeT) == 1){
          pred_y = 0
          for(l in 1 :length(nameMat)){
            pred_y = pred_y + sortie[[j]][1,which(names(sortie[[j]])==nameMat[l])]*seq(X[1],X[length(X)],length.out = min(10,length(X)))^(deg-l+1)
          }
          lines(seq(X[1],X[length(X)],length.out = min(10,length(X))),pred_y, lwd=3.5,lty=3)

        }else{
          abline(v=c(sortie[[j]]$TimeT[-1]),lwd=3.5,lty=3)
          Tt = sortie[[j]]$TimeT[-1]
          for(i in 1 : length(sortie[[j]]$TimeT)){
            if(i == 1){
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[[j]][i,which(names(sortie[[j]])==nameMat[l])]*seq(X[1],Tt[i],length.out = min(10,length(X[which(X <= Tt[i])])))^(deg-l+1)
              }
              lines(seq(X[1],Tt[i],length.out = min(10,length(X[which(X <= Tt[i])]))), pred_y, lwd=3.5,lty=3)
            }else if(i == length(sortie[[j]]$TimeT) ){
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[[j]][i,which(names(sortie[[j]])==nameMat[l])]*seq(Tt[i-1],X[length(X)],length.out = min(10,length(X[which(X >= Tt[i-1])])))^(deg-l+1)
              }
              lines(seq(Tt[i-1],X[length(X)],length.out = min(10,length(X[which(X >= Tt[i-1])]))),pred_y, lwd=3.5,lty=3)
            }else{
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[[j]][i,which(names(sortie[[j]])==nameMat[l])]*seq(Tt[i-1],Tt[i],length.out = min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))])))^(deg-l+1)
              }
              lines(seq(Tt[i-1],Tt[i],length.out = min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))]))),pred_y, lwd=3.5,lty=3)
            }
          }
        }
      }
    }

  }else{
      nameMat = c()
      for(i in 1:(deg+1)){
       if((deg-i+1) == 0){
         nameMat = c(nameMat,"1")
        }else{
         nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
        }
      }

      colnames(list_res[[length(list_res)]]$parameters) = nameMat
      sortie = cbind.data.frame(data.frame("TimeT" = c(X[1],list_res[[length(list_res)]]$TimeTX )),
                                list_res[[length(list_res)]]$parameters,"sigma2" = list_res[[length(list_res)]]$Variances,
                                row.names=as.factor(c(1:(length(list_res[[length(list_res)]]$TimeTX )+1))))

      if(plotG==TRUE){
        plot(X,Y,pch=20,cex.main = 3, cex.lab = 2,cex.axis = 2)
        if(length(sortie$TimeT) == 1){
          pred_y = 0
          for(l in 1 :length(nameMat)){
            pred_y = pred_y + sortie[1,which(names(sortie)==nameMat[l])]*seq(X[1],X[length(X)],length.out = min(10,length(X)))^(deg-l+1)
          }
          lines(seq(X[1],X[length(X)],length.out = min(10,length(X))),pred_y, lwd=3.5,lty=3)

        }else{
          abline(v=c(sortie$TimeT[-1]),lwd=3.5,lty=3)
          Tt = sortie$TimeT[-1]
          for(i in 1 : length(sortie$TimeT)){
            if(i == 1){
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(X[1],Tt[i],length.out =  min(10,length(X[which(X <= Tt[i])])))^(deg-l+1)
              }
              lines(seq(X[1],Tt[i],length.out =  min(10,length(X[which(X <= Tt[i])]))),pred_y, lwd=3.5,lty=3)
            }else if(i == length(sortie$TimeT) ){
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(Tt[i-1],X[length(X)],length.out =  min(10,length(X[which(X >= Tt[i-1])])))^(deg-l+1)
              }
              lines(seq(Tt[i-1],X[length(X)],length.out =  min(10,length(X[which(X >= Tt[i-1])]))),pred_y, lwd=3.5,lty=3)
            }else{
              pred_y = 0
              for(l in 1 :length(nameMat)){
                pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(Tt[i-1],Tt[i],length.out =  min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))])))^(deg-l+1)
              }
              lines(seq(Tt[i-1],Tt[i],length.out =  min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))]))),pred_y, lwd=3.5,lty=3)
            }
          }
        }
      }
    }

  return(sortie)
}



