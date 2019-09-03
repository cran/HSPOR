
#' @title Inference method for two regimes

#' @description H2SPOR is an inference method that estimates, under regularity constraint, the parameters of a
#' polynomial regression model with 2 regimes.
#'
#' @param X A numerical vector corresponding to the explanatory variable. X must be sorted in ascending order
#' if this is not the case, X will be sorted in the function and the corresponding permutation will be applied to Y. The
#' user will be notified by a warning message. In addition, if X contains NAs, they will be deleted from the data and the user will be notified by a warning message.
#' Finally, if X contains duplicate data, the excess data will be deleted and the user will be notified by a warning message.
#' @param Y A numerical vector corresponding to the variable to be explain. It should contain two regimes that
#' could be modelled by polynomials. In addition, if Y contains NAs they will be deleted from the data and the
#' user will be notified by a warning message. Finally, if X contains dupplicate data, the excess data will be deleted and
#' the value of the remaining Y will become the average of the Ys, calculated for this value of X.
#' @param deg The degree of polynomials. The size of X and Y must be greater than 2(deg+2) + 1.
#' @param constraint Number that determines the regularity assumption that is applied for the parameters estimation.
#' By default, the variable is set to 1, i. e. the parameters estimation is done under continuity constraint.
#' If the variable is 0 or 2, the estimation of the parameters will be done without assumption of regularity
#' (constraint = 0) or under assumption of differentiability (constraint = 2). Warning, if the differentiability
#' assumption is not verified by the model, it is preferable not to use it to estimate the model parameters.
#'In addition, if the degree of the polynomials is equal to 1, you cannot use the differentiability assumption.
#' @param EM A Boolean. If EM is TRUE (default), then the function will estimate the parameters
#' of a latent variable polynomial regression model using an EM algorithm. If EM is FALSE then
#' the function will estimate the parameters of the initial polynomial regression model by a fixed point algorithm.
#' @param TimeTrans_Prop A numerical vector. This vector is empty by default. If you want to estimate the model
#' parameters for a fixed jump time value, you can propose this value here.
#' @param plotG A Boolean. If TRUE (default) the estimation results obtained by the H2SPOR function are plotted.
#'
#' @return A dataframe that contains the estimated parameters of the polynomial regression model at two regimes:
#' the jump time, the coefficients of the polynomials and the variances of the two regimes.
#' If plotG = TRUE, the data (X,Y) and the estimated model will be plotted.
#'
#'
#' @import graphics
#' @import corpcor
#' @import stats
#' @export
#'
#'
#' @examples
#'
#'
#' #generated data with two regimes
#' set.seed(1)
#'
#' xgrid1 = seq(0,10,length.out=6)
#' xgrid2 = seq(10.2,20,length.out=6)
#' ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,3)
#' ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,3)
#' xgrid = c(xgrid1,xgrid2)
#' ygrid = c(ygrid1,ygrid2)
#'
#' #Inference of a polynomial regression model with two regimes on these data.
#' #The degree of the polynomials is fixed to 2 and the parameters are estimated
#' #under continuity constraint.
#' H2SPOR(xgrid,ygrid,2,1,EM=FALSE,c())
#'
#' \donttest{
#' set.seed(1)
#' xgrid1 = seq(0,10,by=0.2)
#' xgrid2 = seq(10.2,20,by=0.2)
#' ygrid1 = xgrid1^2-xgrid1+1+ rnorm(length(xgrid1),0,3)
#' ygrid2 = rep(91,length(xgrid2))+ rnorm(length(xgrid2),0,3)
#' xgrid = c(xgrid1,xgrid2)
#' ygrid = c(ygrid1,ygrid2)
#'
#' #Inference of a polynomial regression model with two regimes on these data.
#' #The degree of the polynomials is fixed to 2 and the parameters are estimated
#' #under continuity constraint.
#' H2SPOR(xgrid,ygrid,2,1,EM=FALSE,c())
#' #Executed time : 9.69897 secs (intel core i7 processor)
#' }

H2SPOR = function(X,Y,deg,constraint=1,EM=TRUE,TimeTrans_Prop=c(),plotG=TRUE){

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
    warning("X has been ordered")
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

  if(length(X) < 2*(deg+2) + 1){
    warning(paste("X must contain at least",2*(deg+2)+1,"observations to use this degree. The degree was reduced to", floor((length(X)-5)/2),sep=" "))
    deg = floor((length(X)-5)/2)
  }

  if(length(X) < 7 ){
    stop("X must contain at least 7 observations for a degree equal to 1")
  }

  if(length(TimeTrans_Prop) == 0){
    res = .Trichotomic_algorithm_2_SPOR(X,Y,deg,constraint,EM)
    T_opt <- res[[1]]
    mat_param <- res[[2]]
    sigma2 <- res[[3]]
  }else{
    mat_parami <- .Init_param(X,Y,deg,TimeTrans_Prop[1])[[1]]
    sigma2i <- .Init_param(X,Y,deg,TimeTrans_Prop[1])[[2]]
    res <- .Optimization_2_SPOR(X,Y,deg,TimeTrans_Prop[1],mat_parami,sigma2i,constraint,EM)
    T_opt <- TimeTrans_Prop[1]
    mat_param <- res[[1]]
    sigma2 <- res[[2]]
  }

  sortie <- cbind.data.frame(data.frame("TimeT" = c(X[1],T_opt)),mat_param,"sigma2" = t(sigma2),row.names=c("1","2"))

  if(plotG == TRUE){
    plot(X,Y,pch=20,cex.main = 3, cex.lab = 2,cex.axis = 2)

    nameMat = c()
    for(i in 1:(deg+1)){
      if((deg-i+1) == 0){
        nameMat = c(nameMat,"1")
      }else{
        nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
      }
    }

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
        if((i == 1) & (i != length(sortie$TimeT)) ){
          pred_y = 0
          for(l in 1 : length(nameMat)){
            pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(X[1],Tt[i],length.out = min(10,length(X[which(X <= Tt[i])])))^(deg-l+1)
          }
          lines(seq(X[1],Tt[i],length.out = min(10,length(X[which(X <= Tt[i])]))),pred_y, lwd=3.5,lty=3)
        }else if(i == length(sortie$TimeT) ){
          pred_y = 0
          for(l in 1 :length(nameMat)){
            pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(Tt[i-1],X[length(X)],length.out = min(10,length(X[which(X >= Tt[i-1])])))^(deg-l+1)
          }
          lines(seq(Tt[i-1],X[length(X)],length.out = min(10,length(X[which(X >= Tt[i-1])]))), pred_y, lwd=3.5,lty=3)
        }else{
          pred_y = 0
          for(l in 1 :length(nameMat)){
            pred_y = pred_y + sortie[i,which(names(sortie)==nameMat[l])]*seq(Tt[i-1],Tt[i],length.out = min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))])))^(deg-l+1)
          }
          lines(seq(Tt[i-1],Tt[i],length.out = min(10,length(X[which( (X >= Tt[i-1]) & (X <= Tt[i]))]))), pred_y, lwd=3.5,lty=3)
        }
      }
    }
  }

  return(sortie)
}



