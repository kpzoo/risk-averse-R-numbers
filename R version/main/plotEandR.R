######################################################################
## Plot EpiFilter estimates for effective R and risk averse E 
######################################################################

# Assumptions
# - uses posterior outputs from EpiFilter of R and E
# Inputs - reproduction number estimate (Rmean) and confidence intervals (Rci), 
# estimate (Emean) and confidence intervals (Eci) of risk averse metric, time (tset),
# string for figure (plotname) and folder path to store results (folres)
# Output - .eps plots of best R and E estimates

plotEandR <- function(Rmean, Rci, Emean, Eci, pR1, pE1, pOm1, tset, plotname, folres){
  # Check lengths
  if (length(Emean) != length(Rmean)){
    print(c('[Rmean Emean] lengths', length(Rmean), length(Emean)))
    stop('Inconsistent incidence and reprod. num vectors')
  }else{
    # Plot of R and E estimates on single figure
    pdf(file=paste0(folres, paste0(c(plotname, "_est"), collapse = ""), '.pdf')) 
    par(mfrow=c(1,1))
    
    # Reprod. num estimates and confidence interval
    plot(tset, Rmean, type = 'l', bty = 'l', lwd = 2, col='blue', ylim = c(0, max(Rci[2,]) + 0.2),
         xlab = "t (days)", ylab = "R(t), E(t)")
    polygon(c(tset, rev(tset)), c(Rci[1,], rev(Rci[2,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Rci[3,], rev(Rci[4,])), 
            col =  adjustcolor("dodgerblue", alpha.f = 0.30), border = NA)
    
    # Risk averse estimates and confidence interval
    lines(tset, Emean, type = 'l', bty = 'l', lwd = 2, col='red')
    polygon(c(tset, rev(tset)), c(Eci[1,], rev(Eci[2,])), 
            col =  adjustcolor("indianred", alpha.f = 0.20), border = NA)
    polygon(c(tset, rev(tset)), c(Eci[3,], rev(Eci[4,])), 
            col =  adjustcolor("indianred", alpha.f = 0.30), border = NA)
    lines(tset, rep(1, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    dev.off()
    
    # Single plot of p(R > 1)
    pdf(file=paste0(folres, paste0(c(plotname, "_R1"), collapse = ""), '.pdf')) 
    par(mfrow=c(1,1))
    plot(tset, pR1, type = 'l', bty = 'l', lwd = 2, col='blue', ylim = c(0, 1.05),
         xlab = "t (days)", ylab = "p(R(t) > 1), p(E(t) > 1)")
    lines(tset, pE1, lwd = 2, col = 'red')
    lines(tset, rep(0.5, length(tset)), lwd = 2, col = 'black', lty = 'dashed')
    dev.off()
  }
}