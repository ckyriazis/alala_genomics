library(plyr)

setwd("~/Documents/SDZWA/alala_genomics/analysis/pixy")
data <- read.table("pixy_pi_living.txt", header = T)

mean(data$avg_pi, na.rm=T)


# histogram of pi before filtering
hist(data$avg_pi, xlim = c(0,0.0015), breaks =1000)

# filter out windows with low # of sites
hist(data$no_sites)
data_filtered <- data[data$no_sites>40000,]
plot(data_filtered$no_sites, data_filtered$avg_pi)

# new histogram of pi 
pdf("~/Documents/SDZWA/alala_genomics/analysis/pixy/plots/pi_hist.pdf", width=5, height=4)
hist(data_filtered$avg_pi, xlim = c(0,0.0015), breaks =50, main = "pi in 100kb windows",
     xlab = "pi", col = "chartreuse4")
dev.off()

# filter down to 5% tail of low pi 
quantile(data_filtered$avg_pi, probs = c(0.05, 0.95))
data_filtered_lowpi <- data_filtered[data_filtered$avg_pi<6.056835e-05,]


# how many windows per chrom? is it just driven by chrom length?
chrom_counts <- count(data_filtered_lowpi$chromosome)
colnames(chrom_counts) <- c("Chrom", "freq")

chrom_lengths <- read.csv("chrom_lengths.csv")
chrom_merged <- merge(chrom_counts, chrom_lengths)

pdf("~/Documents/SDZWA/alala_genomics/analysis/pixy/plots/chrom_length_lowPi.pdf", width=5, height=5)
chrom_merged$Length_Mb <- chrom_merged$Length/1e6
plot(chrom_merged$Length_Mb,chrom_merged$freq, pch=19, col = 'chartreuse4', ylab="# of low pi windows", xlab = "Chromosome length (Mb)")
lm <- lm(chrom_merged$freq ~ chrom_merged$Length_Mb)
summary(lm)
abline(a=lm$coefficients[1], b=lm$coefficients[2])
dev.off()


pdf("~/Documents/SDZWA/alala_genomics/analysis/pixy/plots/chrom_pi.pdf", width=6, height=8)
par(mfrow= c(5,3), mar=c(4,2,2,2))

ymax=0.0012 
plot(data$window_pos_1[data$chromosome=="NC_063213.1"], data$avg_pi[data$chromosome=="NC_063213.1"], type = 'l', 
     ylab ="pi", xlab = "position" , col = "chartreuse4", main = "Chromosome 1", ylim=c(0,ymax))
plot(data$window_pos_1[data$chromosome=="NC_063214.1"], data$avg_pi[data$chromosome=="NC_063214.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 2", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063215.1"], data$avg_pi[data$chromosome=="NC_063215.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 3", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063216.1"], data$avg_pi[data$chromosome=="NC_063216.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 4", ylim=c(0,ymax))
plot(data$window_pos_1[data$chromosome=="NC_063217.1"], data$avg_pi[data$chromosome=="NC_063217.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 5", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063218.1"], data$avg_pi[data$chromosome=="NC_063218.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 6", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063219.1"], data$avg_pi[data$chromosome=="NC_063219.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 7", ylim=c(0,ymax))
plot(data$window_pos_1[data$chromosome=="NC_063220.1"], data$avg_pi[data$chromosome=="NC_063220.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 8", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063221.1"], data$avg_pi[data$chromosome=="NC_063221.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 9", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063222.1"], data$avg_pi[data$chromosome=="NC_063222.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 10", ylim=c(0,ymax))
plot(data$window_pos_1[data$chromosome=="NC_063223.1"], data$avg_pi[data$chromosome=="NC_063223.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 11", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063224.1"], data$avg_pi[data$chromosome=="NC_063224.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 12", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063225.1"], data$avg_pi[data$chromosome=="NC_063225.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 13", ylim=c(0,ymax))
plot(data$window_pos_1[data$chromosome=="NC_063226.1"], data$avg_pi[data$chromosome=="NC_063226.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 14", ylim=c(0,ymax), yaxt='n')
plot(data$window_pos_1[data$chromosome=="NC_063227.1"], data$avg_pi[data$chromosome=="NC_063227.1"], type = 'l', 
     ylab ="", xlab = "position" , col = "chartreuse4", main = "Chromosome 15", ylim=c(0,ymax), yaxt='n')

dev.off()



# get mean pi for each chromosome
chroms <- unique(data$chromosome)
chrom_pi <- c()
for(chrom in chroms){
  chrom_pi <- c(chrom_pi, mean(data$avg_pi[data$chromosome==chrom], na.rm = T))
}

# plot mean pi against chromosome length
pdf("~/Documents/SDZWA/alala_genomics/analysis/pixy/plots/chrom_pi_length.pdf", width=4, height=4)
plot(log10(chrom_lengths$Length), chrom_pi, pch=19, col = "chartreuse4", xlab="log10(Chromosome Length)", ylab = "average pi")
lm1 <- lm(chrom_pi~log10(chrom_lengths$Length))
summary(lm1)
abline(a=lm1$coefficients[1], b=lm1$coefficients[2])
dev.off()





