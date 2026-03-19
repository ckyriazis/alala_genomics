

setwd("~/Documents/SDZWA/alala_genomics/simulations/data")

burnin <- read.csv("alala_burnin.csv")
burnin_mu15e8 <- read.csv("alala_burnin_mu15e8.csv")


pdf("~/Documents/SDZWA/alala_genomics/simulations/plots/burnin.pdf", width = 5, height=7)
par(mfrow=c(5,2), mar=c(4,4.5,1,1))
plot(burnin$gen, burnin$popSizep1, type = 'l', col = "chartreuse4", lwd=3, ylab="Population size", xlab = "Generation", ylim=c(50000,90000), xlim=c(0,20000), main = "mu=5e-9")
plot(burnin_mu15e8$gen, burnin_mu15e8$popSizep1, type = 'l', col = "chartreuse4", lwd=3, ylab="Population size", xlab = "Generation", ylim=c(50000,90000), xlim=c(0,20000), main="mu=1.5e-8")


plot(burnin$gen, burnin$FROH_1Mbp1, type = 'l', col = "chartreuse4", lwd=3, ylab=expression(F[ROH]), xlab = "Generation", xlim=c(0,20000), ylim=c(0,0.8))
plot(burnin_mu15e8$gen, burnin_mu15e8$FROH_1Mbp1, type = 'l', col = "chartreuse4", lwd=3, ylab=expression(F[ROH]), xlab = "Generation", xlim=c(0,20000), ylim=c(0,0.8))

plot(burnin$gen, burnin$meanFitnessp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Fitness", xlab = "Generation", ylim=c(0.92,1), xlim=c(0,20000))
plot(burnin_mu15e8$gen, burnin_mu15e8$meanFitnessp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Fitness", xlab = "Generation", ylim=c(0.92,1), xlim=c(0,20000))


plot(burnin$gen, burnin$B_genp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Inbreeding load", xlab = "Generation", ylim=c(0,12), xlim=c(0,20000))
plot(burnin_mu15e8$gen, burnin_mu15e8$B_genp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Inbreeding load", xlab = "Generation", ylim=c(0,12), xlim=c(0,20000))

plot(burnin$gen, burnin$maxLetp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Max lethal freq", xlab = "Generation", ylim=c(0,0.02), xlim=c(0,20000))
plot(burnin_mu15e8$gen, burnin_mu15e8$maxLetp1, type = 'l', col = "chartreuse4", lwd=3, ylab="Max lethal freq", xlab = "Generation", ylim=c(0,0.02), xlim=c(0,20000))



dev.off()