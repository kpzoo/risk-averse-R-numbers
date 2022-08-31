#################################################################################
## Compute effective (R) and risk averse reproduction numbers (E) on empirical data 
## Case study from the Delta COVID-19 wave in Israel. This is from the paper:
## KV Parag and U Obolski, Risk averse reproduction numbers improve resurgence detection,
## medRxiv, 2022:.
## More extensive analyses at https://github.com/kpzoo/risk-averse-R-numbers  
## Easily modified to user-defined data provided local incidence curves are available
#################################################################################

# Functionality shown in this vignette
# - assumes local groups with reproduction numbers Rj and want global statistics
# - load COVID daily incidence curves by testing data from cities in Israel
# - data obtained from dashboard https://datadashboard.health.gov.il/COVID-19/general
# - use serial interval from Nishiura et al 2020 doi:10.1016/j.ijid.2020.02.060
# - estimate R and E using EpiFilter but with sampling from local Rj estimates
# - assumes data is not delayed and any under-reporting is constant
# - uses same serial interval for all groups but this can be modified easily

# Clean the workspace and console
closeAllConnections(); rm(list=ls()); cat("\014")
graphics.off(); start_time = as.numeric(Sys.time())

# Moving average package
library(caTools)

# Set working directory to source
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# Folder path for results
folres = paste0("./results/israel/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

# Load local incidence data by date of testing
Ifiles = dir('data/israel/'); nLoc = length(Ifiles)
# Size of incidence curves - assumes data over same length
Iex = read.csv('data/israel/citycode2610.csv')
tdate = Iex$date; nday = length(tdate)

# Local incidence curve for all groups
Iloc = matrix(0, nLoc, nday)
for (i in 1:nLoc) {
  Iex = read.csv(paste0('data/israel/', Ifiles[i]))
  Iloc[i,] = Iex$new_cases
  # Add weekly smoothing
  Iloc[i,] = round(runmean(Iloc[i,], 7))
}

# Truncate curves so they start from non zero values
ndays = rep(0, nLoc); istarts = ndays; tdays = list()
for (i in 1:nLoc) {
  # Starting id for this group
  iex = which(Iloc[i,] != 0); istarts[i] = iex[1]
  ndays[i] = nday - istarts[i] + 1
  # Days that are valid for estimation
  tdays[[i]] = istarts[i]:nday
}

# Define serial interval from Nishiura et al 2020 (days)
mean_si = 4.7; sd_si = 2.9; pms = c(0,0)
pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si
# Serial interval over max time series length
wdist = pgamma(1:nday, shape = pms[1], scale = pms[2]) - 
  pgamma(0:(nday-1), shape = pms[1], scale = pms[2])

# Total infectiousness of each local group
Lloc = matrix(0, nLoc, nday); Lloc[,1] = Iloc[,1]
for(j in 1:nLoc){
  for(i in 2:nday){
    # Total infectiousness
    Lloc[j,i] = sum(Iloc[j, seq(i-1, 1, -1)]*wdist[1:(i-1)])    
  }
}

# Sum of infections and total infections 
Itot = colSums(Iloc); Ltot = colSums(Lloc)

######################################################################
## Local scale reproduction number Rj estimates
######################################################################

# Setup grid [Rmin Rmax] noise eta 
Rmin = 0.01; Rmax = 10; eta = 0.1
# Uniform prior over grid of size m
m = 1000; pR0 = (1/m)*rep(1, m);conf = 0.025

# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)
# Store outputs of mean Rj, 95% CIs and P(Rj > 1)
Rmeanloc = list(); Rciloc = list(); pRloc1 = list()

# Sample from local Rj posteriors
nSamp = 10000; Rjsamp = vector("list", nLoc)

# Estimate local reproduction numbers with EpiFilter
for (j in 1:nLoc) {
  # Filtered estimates as list [Rmed, Rci, Rmean, pR, pRup, pstate]
  Rfilt = epiFilter(Rgrid, m, eta, pR0, ndays[j], Lloc[j,tdays[[j]]], Iloc[j,tdays[[j]]], conf)
  # Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
  Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, ndays[j], Rfilt$pstate, conf)
  
  # Sample from posterior distribution at each time
  Rsamp = matrix(0, ndays[j], nSamp)
  for (i in 1:ndays[j]) {
    Rsamp[i,] = sample(Rgrid, nSamp, replace = TRUE, prob = Rsmooth$qR[i,])
  }
  # Store samples across time from each location
  Rjsamp[[j]] = Rsamp
  # Save estimates for constructing consensus statistics
  Rmeanloc[[j]] = Rsmooth$Rmean; Rciloc[[j]] = Rsmooth$Rci 
}

######################################################################
## Compute consensus statistics E and standard R
######################################################################

# Estimate global effective reproduction number
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Ltot[1:nday], Itot[1:nday], conf)
Rsmooth = epiSmoother(Rgrid, m, Rfilt$pR, Rfilt$pRup, nday, Rfilt$pstate, conf)
R = Rsmooth$Rmean; Rci = Rsmooth$Rci 
pR1 = rep(0, nday); id1 = which(Rgrid >= 1); id1 = id1[1]
for (i in 1:nday){
  pR1[i] = sum(Rsmooth$qR[i, id1:m])
}

# Save effective global R estimates
nam1 = paste0(c(folres, "Rmean"), collapse = '')
nam2 = paste0(c(folres, "Rci"), collapse = '')
nam3 = paste0(c(folres, "pR1"), collapse = '')
write.table(R, file = nam1, row.names=FALSE, col.names=FALSE)
write.table(Rci, file = nam2, row.names=FALSE, col.names=FALSE)
write.table(pR1, file = nam3, row.names=FALSE, col.names=FALSE)

# Estimate risk averse reproduction number E from Rj samples
Esamp = matrix(0, nday, nSamp)
for (j in 1:nday) {
  # Collect samples across groups for a given day
  Rjday = list(); id = 1
  for (i in 1:nLoc){
    # Local samples from every Rj posterior
    if(j >= istarts[i]){
      Rjday[[id]] = Rjsamp[[i]][j - istarts[i]+1,]
      id = id + 1
    }
  }
  # E is sum of squares divided by sum of samples
  nLocActive = length(Rjday); Ej = matrix(0, nLocActive, nSamp)
  for(i in 1:nLocActive){
    # Individual contributions to E from groups
    Ej[i,] = Rjday[[i]]
  }
  # Consensus statistic E
  num = colSums(Ej^2); den = colSums(Ej)
  Esamp[j,] = num/den
}

# Mean and CIs of risk averse reproduction number
Emean = rowMeans(Esamp); Eci = matrix(0, 4, nday)
for (i in 1:nday) {
  Eci[, i] = quantile(Esamp[i,], c(0.025, 0.975, 0.25, 0.75))
}
# Prob of E > 1 (resurgence)
pE1 = rep(0, nday)
for (i in 1:nday){
  pE1[i] = length(which(Esamp[i,] > 1))/nSamp
}

# Save E estimates
nam1 = paste0(c(folres, "Emean"), collapse = '')
nam2 = paste0(c(folres, "Eci"), collapse = '')
nam3 = paste0(c(folres, "pE1"), collapse = '')
write.table(Emean, file = nam1, row.names=FALSE, col.names=FALSE)
write.table(Eci, file = nam2, row.names=FALSE, col.names=FALSE)
write.table(pE1, file = nam3, row.names=FALSE, col.names=FALSE)

######################################################################
## Visualisation and output of results (as csv files)
######################################################################


# Plot all estimates against time
plotEandR(R, Rci, Emean, Eci, pR1, pE1, pOm1, 1:nday, 'ER', folres)

# Log running time in minutes
end_time = as.numeric(Sys.time()); exec_time = end_time - start_time 
exec_time = exec_time/60; exec_time = round(exec_time, 3)
print(paste0(c("Completed in ", exec_time, " mins"), collapse = ""))
