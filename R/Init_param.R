
.Init_param <- function(X,Y,deg,timeT){

  nbT <- length(timeT)

  nameMat = c()
  for(i in 1:(deg+1)){
    if((deg-i+1) == 0){
      nameMat = c(nameMat,"1")
    }else{
      nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
    }
  }

  param_init <- matrix(0,nbT+1,(deg+1),dimnames = list(c(),nameMat))
  sigma2 <- matrix(0,1,nbT+1)
  m <- 1

  #We initiate the parameters'models for each regime
  while(m <= (nbT+1)){
    if(m == 1){
      x <- X[which(X <= timeT[m])]
      y <- Y[which(X <= timeT[m])]
    }
    if(m > 1 & m < (nbT+1)){
      x <- X[which(X > timeT[m-1] & X <= timeT[m])]
      y <- Y[which(X > timeT[m-1] & X <= timeT[m])]
    }
    if(m == (nbT+1)){
      x <- X[which(X > timeT[m-1])]
      y <- Y[which(X > timeT[m-1])]
    }

    #Regression
    obs_table <- data.frame()
    for(i in 1:deg){
      obs_table <- as.data.frame(t(rbind.data.frame(t(obs_table),x^i)))
    }

    lmR <- lm(y ~ ., data = obs_table)
    param_init[m,] <- rev(lmR$coefficients)
    sigma2[,m] <- sum((lmR$residuals)^2)/(length(y)-(deg+1))
    m <- m + 1
  }

  list(param_init,sigma2)
}

