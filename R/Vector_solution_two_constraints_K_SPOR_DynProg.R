
.Vector_solution_two_constraints_K_SPOR_DynProg <- function(X,Y,deg,constraint,begin_point,end_point){

  # constraint == 1 (continuity)
  if(constraint == 1){
    VsolutionC <- c(begin_point[2],end_point[2])
  }

  # constraint == 2 (continuity & derivability)
  if(constraint == 2){
    VsolutionC <- c(begin_point[2],end_point[2],begin_point[3],end_point[3])
  }

  VsolutionM <- rep(0,deg+1)
  for(i in 1 : (deg+1)){
    VsolutionM[i] <- sum((X^(deg+1-i))*Y)
  }

  if(constraint == 0){
    Vsolution <- VsolutionM
  }else{
    Vsolution <- c(VsolutionC,VsolutionM)
  }
  return(Vsolution)
}
