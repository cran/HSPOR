
.Transform_vector_to_matrix <- function(vparam,deg,timeT){

  nbT <- length(timeT)

  nameMat = c()
  for(i in 1:(deg+1)){
    if((deg-i+1) == 0){
      nameMat = c(nameMat,"1")
    }else{
      nameMat = c(nameMat,paste("X^",(deg-i+1),sep=""))
    }
  }

  mat_param <- matrix(vparam[1:((nbT+1)*(deg+1))],nrow <- nbT+1,ncol <- deg+1,dimnames = list(c(),nameMat))

  return(mat_param)
}
