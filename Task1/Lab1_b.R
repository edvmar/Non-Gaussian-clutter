if(!is.null(dev.list())) dev.off()
#rm(list=ls())
#cat("\014")
#setwd("~/OneDrive/Dokument/Skolan/Kurser/Aktiva/ProjectAtSaab/git_code/Non-Gaussian-clutter")
library(kdist)
library(cmvnorm)

#-------Task b)-------

# Here we have:
#   Gaussian detector
#   Compound Gaussian clutter

n = 1000 # Number of simulations in Monte Carlo
alpha = 2
theta = 0
s <- complex(modulus = alpha, argument = theta) # signal 

etas = c(0.1, 0.5, 1, 2, 5, 10) # change later
P_FA = rep(0, length(etas)) # store for plotting 
P_TD = rep(0, length(etas))

for (iEta in 1:length(etas)){
  
  eta = etas[iEta]
  
  # Why new sample for each eta in description? 
  a = rk(n, shape = 1, scale = 1, intensity = FALSE)
  b = rk(n, shape = 1, scale = 1, intensity = FALSE)
  c = complex(real = a, imaginary = b)
  
  r_0 = c
  r_1 = c + s
  
  # TODO: Find these correctly from PDF for CN 
  # TODO: Move to for loop? Don't need to store them.. 
  f0_r0 = rep(0,n)
  f1_r0 = rep(0,n)
  
  f0_r1 = rep(0,n)
  f1_r1 = rep(0,n)
  
  #for(i in 1:n){
  #  f_0[i] <- dcmvnorm(r_0[i], 0, 1)
  #  f_0[i] <- dcmvnorm(r_1[i], 0, 1)
  #}
  
  for (j in 1:n){
    
    #False Alarm 
    if(f1_r0[j]/f0_r0[j] > eta){ # False detection 
      P_FA[iEta] = P_FA[iEta] + f0_r0[j]
    }
    #True Detection
    if(f1_r1[j]/f0_r1[j] > eta){ # True detection 
      P_TD[iEta] = P_TD[iEta] + f1_r1[j]
    }
  }

}

P_FA = P_FA/N.  #Normalize 
P_TD = P_TD/N


# TODO: Plot P_FA vs P_TD correctly with log scales etc
P_FA = c(0, 0.2, 0.3, 0.4)
P_TD = c(0, 0.4, 0.6, 0.8)
plot(P_FA, P_TD, type = "b", frame = FALSE, pch = 19, 
     col = "black", xlab = "P_FA", ylab = "P_TD")
legend("right", legend=c("SIR = ? "),
       col=c("black"), lty = 1:2, cex=0.8)












