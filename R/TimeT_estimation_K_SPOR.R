
.TimeT_estimation_K_SPOR <- function(X,Y,deg,grid,current_K,K,constraint,EM,i=1,timeT=c(),mat_res=list(),j_val=c()){

  resFi <- c()

  if(current_K != K){
    if(current_K == 1){
      beginning <- 1
    }else{
      beginning <-  (i+(deg+2))
    }

    mat_res[[current_K]] <- matrix(0,nrow=(K-current_K),ncol=length(grid[[current_K]]))

    for(j in (beginning):(length(grid[[current_K]]))){
      j_val[current_K] = j

      timeT[current_K] <- grid[[current_K]][j]

      if(current_K == (K-1)){
        res <-  .Optimization_K_SPOR(X,Y,deg,timeT,EM,constraint)
        mat_res[[current_K]][1,j] <- res[[3]]
      }else{
        res <-  .TimeT_estimation_K_SPOR(X,Y,deg,grid,current_K+1,K,constraint,EM,j,timeT,mat_res,j_val)
        mat_res <- res$mat_res
      }
    }

    if(current_K > 1){
      for(k in 1:(nrow(mat_res[[current_K-1]])-1)){
        mini <- which.min(mat_res[[current_K]][1,beginning:length(grid[[current_K]])])
        mat_res[[current_K-1]][k,j_val[current_K-1]] <- mat_res[[current_K]][k,beginning+mini-1]
      }
      mat_res[[current_K-1]][nrow(mat_res[[current_K-1]]),j_val[current_K-1]] <- grid[[current_K]][beginning+mini-1]
    }
    if(current_K == 1){
      for(l in 1:nrow(mat_res[[current_K]])){
        miniFi <- which.min(mat_res[[current_K]][1,])
        resFi <- c(resFi,mat_res[[current_K]][l,miniFi])
      }
      resFi <- c(resFi, grid[[current_K]][miniFi])
    }
  }
  list("mat_res" = mat_res, "resFi" = resFi)
}
