
.TimeT_grid_K_SPOR <- function(X,deg,K){
  grid_Tt <- list()
  for( i in 1 : (K-1)){
    grid_Tt[[i]] <- sapply((deg+2):(length(X)- ((deg+2)*(K-i))), function(j) (X[j+1]+X[j])/2)
  }
  return(grid_Tt)
}
