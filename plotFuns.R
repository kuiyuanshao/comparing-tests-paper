#######################################################################################################################
############################################### Plotting Functions ####################################################
#######################################################################################################################

plotfun <- function(exp_pmatrix, pmatrix, exp_pmatrix2, pmatrix2){
  for (i in 1:5){
    par(pty="s")
    mym <- matrix(1:2, ncol = 2)
    layout(mym)
    qqplot(x = 0, y = 0, 
           xlim = c(0, 
                    max(-log10(exp_pmatrix[i, ]), 
                        -log10(exp_pmatrix[i + 5, ]))),
           ylim = c(0, 3), 
           col = rgb(red = 1, green = 1, blue = 1, alpha = 0), 
           xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
    log_p <- -log10(pmatrix[i, ])
    log_p2 <- -log10(pmatrix[i + 5, ])
    log_p[which(log_p > 2.5)] <- 2.5
    log_p2[which(log_p2 > 2.5)] <- 2.5
    points(x = -log10(exp_pmatrix[i, ]), 
           y = log_p, col = "#4896E0")
    points(x = -log10(exp_pmatrix[i + 5, ]), 
           y = log_p2, col = "#CF5D6D")
    legend("topleft", col = c("#4896E0", "#CF5D6D"), 
           c("F-based", "Chi-square based"), 
           pch = 1, bty = "n")
    abline(0, 1)
    
    
    qqplot(x = 0, y = 0, 
           xlim = c(0, 
                    max(-log10(exp_pmatrix2[i, ]), 
                        -log10(exp_pmatrix2[i + 5, ]))),
           ylim = c(0, 3), 
           col = rgb(red = 1, green = 1, blue = 1, alpha = 0), 
           xlab = "Expected -log10(p)", ylab = "Observed -log10(p)")
    log_p <- -log10(pmatrix2[i, ])
    log_p2 <- -log10(pmatrix2[i + 5, ])
    log_p[which(log_p > 2.5)] <- 2.5
    log_p2[which(log_p2 > 2.5)] <- 2.5
    points(x = -log10(exp_pmatrix2[i, ]), 
           y = log_p, col = "#4896E0")
    points(x = -log10(exp_pmatrix2[i + 5, ]), 
           y = log_p2, col = "#CF5D6D")
    abline(0, 1)
  }
}