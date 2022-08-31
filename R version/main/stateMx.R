######################################################################
## State model choices for EpiFilter 2.0
######################################################################

# Choices - reproduction number state space model can be selected from various choices
# Inputs - reproduction number grid (Rgrid), size of grid (m), state type with parameters (stateMod)
# Output - state transition matrix (stateMx)

stateMx <- function(Rgrid, m, stateMod){

  # Identification of specific model
  stateID = stateMod$id
  # Instantiate state matrix
  pstate = matrix(0, m, m)
  
  # Compute state matrix based on chosen prior model
  if(stateID == 0){
    # Two scaling parameters (first is eta)
    param = stateMod$params
    
    # Standard EpiFilter (scaled diffusion) with option to remove scaling
    for(j in 1:m){
      pstate[j, ] = dnorm(Rgrid[j], Rgrid, sqrt(Rgrid)*param[1] + param[2])
    }
  }
  if(stateID == 1){
    # One scaling Rayleigh parameter
    param = stateMod$params
    
    # Rayleigh constraint on state
    for(j in 1:m){
      #pstate[j, ] = drayleigh(Rgrid[j], scale = sqrt(Rgrid)*param[1], log = FALSE)
      pstate[j, ] = drayleigh(Rgrid[j], scale = param[1], log = FALSE)
    }
  }
  if(stateID == 2){
    
  }
  
  stateMx = t(pstate)
}