if(!is.null(dev.list())) dev.off()
rm(list=ls())
cat("\014")
setwd("~/OneDrive/Dokument/Skolan/Kurser/Aktiva/ProjectAtSaab/git_code/Non-Gaussian-clutter")
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
eta = 3
a = rk(n, shape = 1, scale = 1, intensity = FALSE)
b = rk(n, shape = 1, scale = 1, intensity = FALSE)
c = complex(real = a, imaginary = b)

r_0 = rep(0, n)

for(i in 1:n){
  r_0[i] <- rcmvnorm(1, c[i], 1)
}



y0 = c






  
















#-------Task c)-------

# Here we have:
#   Compound Gaussian detector
#   Gaussian clutter


#-------Task d)-------

# Here we have:
#   Compound Gaussian detector
#   Compound Gaussian clutter