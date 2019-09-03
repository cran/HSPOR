
.H2SPOR_DynProg_rec <- function(X,Y,deg,nb_saut,n_left,n_right,LineToChange,V_Tt,matP,var,dir,constraint,pt_const,BIC_ok,EM){

  ecart <- 1
  pct_ecart = 13

  if((length(X) > ((deg+2)*2)+1) & (ecart > 0) & (pct_ecart>=10)){

    if(dir == 'left'){
      n_left = n_left+1
    }else{
      n_right = n_right+1
    }

    if(nb_saut==1){
      res = .Trichotomic_algorithm_one_constraint_2_SPOR_DynProg(X,Y,deg,constraint,EM,0.001)
    }else{
      res = .Trichotomic_algorithm_many_constraints_2_SPOR_DynProg(X,Y,deg,constraint,EM,pt_const,dir,0.001)
    }

    T_opt <- res[[1]]
    #print(T_opt)
    mat_param <- res[[2]]
    sigma2 <- res[[3]]
    MLL <- res[[4]]

    #calcul du BIC pour valider le saut
    if(constraint == 1){
      nb_param <- (2*(deg+1)+1) - (nrow(pt_const)+1)
    }else if(constraint == 2){
      nb_param <- (2*(deg+1)+1) - 2*(nrow(pt_const)+1)
    }else{
      nb_param <- 2*(deg+1)+1
    }

    BIC_saut <- 2*MLL + log(length(X))*(nb_param)

    #BIC without jump
    obs_table <- data.frame()
    for(i in 1:deg){
      obs_table <- as.data.frame(t(rbind.data.frame(t(obs_table),X^i)))
    }
    lmR <- lm(Y ~ ., data = obs_table)
    par <- rev(lmR$coefficients)
    sig2 <- sum((lmR$residuals)^2)/(length(Y)-(deg+1))
    s = 0
    for(i in 1:(deg+1)){
      s = s + (par[i]*(X^(deg+1-i)))
    }

    dens<- (1/sqrt(2*pi*sig2)) * exp( - ((Y - s)^2)/(2*sig2))

    BIC_ss_saut <- -2*sum(log(dens)) + log(length(X))*(deg+1)

    ecart = BIC_ss_saut-BIC_saut
    #print(paste("ecart",ecart,sep=" "))

    pct_ecart = ecart/BIC_ss_saut*100
    #print(paste("pct_ecart",pct_ecart,sep=" "))
    if((ecart > 0) & (pct_ecart>=10)){

      V_Tt <- c(V_Tt,T_opt)

      if(nb_saut ==1){
        BIC_ok = matrix(c(BIC_saut,BIC_ss_saut),nrow=1,ncol=2)
      }else{
        BIC_ok = rbind(BIC_ok,matrix(c(BIC_saut,BIC_ss_saut),nrow=1,ncol=2))
      }
      BIC_ok = BIC_ok[order(V_Tt),]
      #print(paste("BIC_ok",BIC_ok,sep=" "))

      V_Tt <- V_Tt[order(V_Tt)]
      #print(paste("V_T",V_Tt,sep=" "))

      if(nb_saut == 1){
        matP <- mat_param
        var <- sigma2
      }else{
        LineToChange <- which(V_Tt == T_opt)
        #print(nrow(matP))
        #print(mat_param)
        if(LineToChange == 1){
          matP <- rbind(mat_param,matP[2:nrow(matP),])
          var <- matrix(c(sigma2,var[1,(2:ncol(var))]),nrow = 1,ncol=(ncol(var)+1))
        }else if(LineToChange == nrow(matP)){
          matP <- rbind(matP[1:(nrow(matP)-1),],mat_param)
          var <- matrix(c(var[1,1:(ncol(var)-1)],sigma2),nrow=1,ncol=(ncol(var)+1))
        }else{
          matP <- rbind(matP[1:(LineToChange-1),],mat_param,matP[(LineToChange+1):nrow(matP),])
          var <- matrix(c(var[1,1:(LineToChange-1)],sigma2,var[1,(LineToChange+1):ncol(var)]),nrow=1,ncol=(ncol(var)+1))
        }
      }
      #print(matP)

      nb_saut <- nb_saut + 1

      X1 <- X[which(X <= T_opt)]
      Y1 <- Y[which(X <= T_opt)]

      #Calculation of the continuity or differentiability points
      if(nrow(matP) == 2){
        sum_cont <- 0
        sum_der <- 0
        for(i in 1:(deg+1)){
          sum_cont <- sum_cont + matP[2,i]*T_opt^(deg+1-i)
          sum_der <- sum_der + ((deg+1-i)*matP[2,i])*T_opt^(deg-i)
        }
        pt_const <- matrix(c(T_opt,sum_cont,sum_der),nrow=1,ncol=3)
      }else{
        if(which(V_Tt == T_opt) == 1){
          sum_cont <- 0
          sum_der <- 0
          for(i in 1:(deg+1)){
            sum_cont <- sum_cont + matP[2,i]*T_opt^(deg+1-i)
            sum_der <- sum_der + ((deg+1-i)*matP[2,i])*T_opt^(deg-i)
          }
          pt_const <- matrix(c(T_opt,sum_cont,sum_der),nrow=1,ncol=3)
        }else{
          sum_cont1 <- 0
          sum_der1 <- 0
          sum_cont2 <- 0
          sum_der2 <- 0

          for(i in 1:(deg+1)){
            sum_cont1 <- sum_cont1 + matP[(which(V_Tt == T_opt)-1),i]*V_Tt[(which(V_Tt == T_opt)-1)]^(deg+1-i)
            sum_der1 <- sum_der1 + ((deg+1-i)*matP[(which(V_Tt == T_opt)-1),i])*V_Tt[(which(V_Tt == T_opt)-1)]^(deg-i)
            sum_cont2 <- sum_cont2 + matP[(which(V_Tt == T_opt)+1),i]*T_opt^(deg+1-i)
            sum_der2 <- sum_der2 + ((deg+1-i)*matP[(which(V_Tt == T_opt)+1),i])*T_opt^(deg-i)
          }
          pt_const <- matrix(c(V_Tt[(which(V_Tt == T_opt)-1)],sum_cont1,sum_der1,T_opt,sum_cont2,sum_der2),nrow=2,ncol=3,byrow=T)
        }
      }

      resG <- .H2SPOR_DynProg_rec(X1,Y1,deg,nb_saut,n_left,n_right,LineToChange,V_Tt,matP,var,dir="left",constraint,pt_const,BIC_ok,EM)
      V_Tt <- resG$V_Tt
      matP <- resG$matP
      var <- resG$var
      n_left <- resG$n_left
      n_right <- resG$n_right
      LineToChange <- resG$LineToChange
      BIC_ok = resG$BIC_ok

      X2 <- X[which(X > T_opt)]
      Y2 <- Y[which(X > T_opt)]

      if(nrow(matP) == 2){
        sum_cont_r <- 0
        sum_der_r <- 0
        for(i in 1:(deg+1)){
          sum_cont_r <- sum_cont_r + matP[1,i]*T_opt^(deg+1-i)
          sum_der_r <- sum_der_r + ((deg+1-i)*matP[1,i])*T_opt^(deg-i)
        }
        pt_const <- matrix(c(T_opt,sum_cont_r,sum_der_r),nrow=1,ncol=3)
      }else{
        if(which(V_Tt == T_opt) == length(V_Tt)){
          sum_cont_r <- 0
          sum_der_r <- 0
          for(i in 1:(deg+1)){
            sum_cont_r <- sum_cont_r + matP[(length(V_Tt)),i]*T_opt^(deg+1-i)
            sum_der_r <- sum_der_r + ((deg+1-i)*matP[(length(V_Tt)),i])*T_opt^(deg-i)
          }
          pt_const <- matrix(c(T_opt,sum_cont_r,sum_der_r),nrow=1,ncol=3)
        }else{
          sum_cont_r1 <- 0
          sum_der_r1 <- 0
          sum_cont_r2 <- 0
          sum_der_r2 <- 0
          for(i in 1:(deg+1)){
            sum_cont_r1 <- sum_cont_r1 + matP[(which(V_Tt == T_opt)),i]*T_opt^(deg+1-i)
            sum_der_r1 <- sum_der_r1 + ((deg+1-i)*matP[(which(V_Tt == T_opt)),i])*T_opt^(deg-i)
            sum_cont_r2 <- sum_cont_r2 + matP[(which(V_Tt == T_opt)+2),i]*V_Tt[which(V_Tt == T_opt)+1]^(deg+1-i)
            sum_der_r2 <- sum_der_r2 + ((deg+1-i)*matP[(which(V_Tt == T_opt)+2),i])*V_Tt[which(V_Tt == T_opt)+1]^(deg-i)
          }
          pt_const <- matrix(c(T_opt,sum_cont_r1,sum_der_r1,V_Tt[which(V_Tt == T_opt)+1],sum_cont_r2,sum_der_r2),nrow=2,ncol=3,byrow=T)
        }
      }

      resD <- .H2SPOR_DynProg_rec(X2,Y2,deg,nb_saut,n_left,n_right,LineToChange,V_Tt,matP,var,dir="right",constraint,pt_const,BIC_ok,EM)
      V_Tt <- resD$V_Tt
      matP <- resD$matP
      var <- resD$var
      n_left <- resD$n_left
      n_right <- resD$n_right
      LineToChange <- resD$LineToChange
      BIC_ok = resD$BIC_ok

    }else{
      ecart <- 0
    }
  }
  #print(matP)
  list("V_Tt" = V_Tt, "matP" = matP,"var" = var,"n_left" = n_left,"n_right" = n_right,"LineToChange" = LineToChange,"BIC_ok"=BIC_ok)
}
