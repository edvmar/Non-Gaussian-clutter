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

n = 1000  # Number of simulations in Monte Carlo / sample size 
alpha = 2
theta = 0
s <- complex(modulus = alpha, argument = theta) # signal 

etas = c(0.5, 1, 2, 5, 10) # change later
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
  
  
  # TODO: Find these correctly from PDF for CN => I think I did it right now?
  # f( 2+i | CN(0,1) ] = dcmvnorm(c(2,1), c(0,0), diag(2)/2)
 
  for (j in 1:n){
    #False Alarm 
    f0_r0 <- dcmvnorm(c(Re(r_0[j]), Im(r_0[j])), c(0,0), diag(2)/2)             # f(r|H0), r \in r_0 
    f1_r0 <- dcmvnorm(c(Re(r_0[j]), Im(r_0[j])), c(Re(s),Im(s)), diag(2)/2)     # f(r|H1), r \in r_0 

    if(f1_r0/f0_r0 > eta){ # False detection 
      P_FA[iEta] = P_FA[iEta] + f0_r0
    }
    
    #True Detection
    f0_r1 <- dcmvnorm(c(Re(r_1[j]), Im(r_1[j])), c(0,0), diag(2)/2)             # f(r|H0), r \in r_1
    f1_r1 <- dcmvnorm(c(Re(r_1[j]), Im(r_1[j])), c(Re(s),Im(s)), diag(2)/2)     # f(r|H1), r \in r_1 
    
    if(f1_r1/f0_r1 > eta){ # True detection 
      P_TD[iEta] = P_TD[iEta] + f1_r1
    }
  }

}

P_FA = P_FA/n  #Normalize 
P_TD = P_TD/n


# TODO: Plot P_FA vs P_TD correctly with log scales etc
# TODO: Guess we should do this for several SIR:s as well in the same plot..
plot(P_FA, P_TD, type = "b", frame = FALSE, pch = 19, 
     col = "black", xlab = "P_FA", ylab = "P_TD")
legend("right", legend=c("SIR = ? "),
       col=c("black"), lty = 1:2, cex=0.8)













