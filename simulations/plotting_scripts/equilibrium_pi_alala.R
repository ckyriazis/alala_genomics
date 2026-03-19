

setwd("~/Documents/SDZWA/alala_genomics/simulations/data")

burnin <- read.csv("neutral_burnin.csv")

plot(burnin$gen, burnin$meanHetp1)


pi <- mean(tail(burnin$meanHetp1, n = 15))

# Ne at equlibrium when K=1000
Ne <- pi/(4*5e-9*2.31/3.31)
Ne


ratio <- Ne/1000
ratio

33500/ratio


