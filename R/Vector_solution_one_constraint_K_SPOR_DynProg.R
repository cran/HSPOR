
.Vector_solution_one_constraint_K_SPOR_DynProg <- function(X,Y,deg,constraint,constraint_point){

  # constraint == 1 (continuity)
  if(constraint == 1){
    VsolutionC <- constraint_point[2]
  }

  # constraint == 2 (continuity & derivability)
  if(constraint == 2){
    VsolutionC <- c(constraint_point[2],constraint_point[3])
  }

  d <- deg
  VsolutionM <- rep(0,deg+1)
  for(i in 1 : (deg+1)){
    VsolutionM[i] <- sum((X^d)*Y)
    d <- d-1
  }

  if(constraint == 0){
    Vsolution <- VsolutionM
  }else{
    Vsolution <- c(VsolutionC,VsolutionM)
  }

  return(Vsolution)
}
