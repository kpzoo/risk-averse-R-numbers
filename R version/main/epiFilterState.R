######################################################################
## Bayesian recursive filtering with different state models for EpiFilter
# From: Parag, KV (2021). Improved real-time estimation of reproduction numbers at
# low case incidence and between epidemic waves. PLOS Comput Biol 17(9):e1009347.
######################################################################

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model can be input from stateMx.R

# Inputs - reproduction number grid (Rgrid), size of grid (m), diffusion noise (eta), state matrix (pstate)
# prior on R (pR0), max time (nday), total infectiousness (Lday), incidence (Iday), confidence (a)

# Output - mean (Rmean), median (Rmed), 50% and 95% quantiles of estimates (Rhat),
# causal posterior over R (pR), pre-update (pRup) and state transition matrix (pstate)

# Sources of errors 
# - if Rgrid is poorly specified (m too small or Rmin/Rmax too tight) due to Rcdf not being computable
# - if Iday[1] = 0 as the epidemic must be seeded with at least one case (this is a renewal model trait)

epiFilterState <- function(Rgrid, m, eta, pR0, nday, Lday, Iday, a, pstate){
  
  # Probability vector for R and prior
  pR = matrix(0, nday, m); pRup = pR
  pR[1, ] = pR0; pRup[1, ] = pR0
  
  # Mean and median estimates
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rci = matrix(0, 4, nday)
  
  # Initialise mean and CDF of prior 
  Rmean[1] = pR[1, ]%*%Rgrid; Rcdf0 = cumsum(pR0)
  
  # Initialise quartiles
  idm = which(Rcdf0 >= 0.5, 1); Rmed[1] = Rgrid[idm[1]]
  id1 = which(Rcdf0 >= a, 1); id2 = which(Rcdf0 >= 1-a, 1)
  id3 = which(Rcdf0 >= 0.25, 1); id4 = which(Rcdf0 >= 0.75, 1)
  Rci[1, 1] = Rgrid[id1[1]]; Rci[2, 1] = Rgrid[id2[1]]
  Rci[3, 1] = Rgrid[id3[1]]; Rci[4, 1] = Rgrid[id4[1]]
  
  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    rate = Lday[i]*Rgrid
    # Probabilities of observations
    pI = dpois(Iday[i], rate)
    
    # State predictions for R
    pRup[i, ]  = pR[i-1, ]%*%pstate
    # Update to posterior over R
    pR[i, ] = pRup[i, ]*pI
    pR[i, ] = pR[i, ]/sum(pR[i, ])
    
    # Posterior mean and CDF
    Rmean[i] = pR[i, ]%*%Rgrid
    Rcdf = cumsum(pR[i, ])
    
    # Quantiles for estimates
    idm = which(Rcdf >= 0.5, 1); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rci[1, i] = Rgrid[id1[1]]; Rci[2, i] = Rgrid[id2[1]]
    Rci[3, i] = Rgrid[id3[1]]; Rci[4, i] = Rgrid[id4[1]]
  }
  
  # Main outputs: estimates of R and states
  out = list(Rmed, Rci, Rmean, pR, pRup, pstate)
  names(out) = c('Rmed', 'Rci', 'Rmean', 'pR', 'pRup', 'pstate')
  epiFilterState = out
}