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
s <- complex(modulus = alpha, argument = theta)
eta = 0.01
a = rk(n, shape = 1, scale = 1, intensity = FALSE)
b = rk(n, shape = 1, scale = 1, intensity = FALSE)
c = complex(real = a, imaginary = b)

r_0 = c
r_1 = c + s

f_0 = rep(0, n)
f_1 = rep(0, n)

#for(i in 1:n){
#  f_0[i] <- dcmvnorm(r_0[i], 0, 1)
#  f_0[i] <- dcmvnorm(r_1[i], 0, 1)
#}

f_0 = 1/pi*exp(-abs(r_0)^2)
f_1 = 1/pi*exp(-abs(r_1)^2)

# threshhold = 2 # ? 
# P_FA = 0
# P_TD = 0
# for(i in 1:n){
#   if (abs(r_0[i]) > threshhold){
#     P_FA <- P_FA + f_0[i]
#     }
#   if (abs(r_1[i]) > threshhold){
#     P_TD <- P_TD + f_0[i]
#     }
# }
# P_FA = P_FA/n
# P_TD = P_TD/n

P_FA = 0
P_TD = 0
for(i in 1:n){
  LRT = f_1[i]/f_0[i]
  if (LRT > eta){
    P_FA <- P_FA + f_0[i]
    P_TD <- P_TD + f_1[i]
  }
}
P_FA = P_FA/n
P_TD = P_TD/n


























#-------Task c)-------

# Here we have:
#   Compound Gaussian detector
#   Gaussian clutter


#-------Task d)-------

# Here we have:
#   Compound Gaussian detector
#   Compound Gaussian clutter