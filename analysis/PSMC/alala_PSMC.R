


pdf("~/Documents/SDZWA/alala_genomics/analysis/PSMC/plots/alala.pdf", width=7, height=6)
par(mar=c(5,5,2,2))


setwd("~/Documents/SDZWA/alala_genomics/analysis/PSMC/data/175AL_FilterB_Round1_autosomes_PASS")

## read in alala files
alala_files <- list.files(pattern=".txt")
alala_data <- data.frame(matrix(nrow = 0, ncol = 5))

for(rep in seq_along(alala_files)){
  data <- read.table(alala_files[rep], sep="\t", header=F)
  alala_data <- rbind(alala_data,cbind(data,rep))
}

## plot main alala result
rep1 <- alala_data[alala_data$rep == 1, ]
start_year <- 500
lwd <- 10

#full_results_alala <- read.table("~/Documents/SDZWA/alala_genomics/analysis/PSMC/data/175AL_FilterB_Round1_autosomes_PASS_AL95.psmc.plot.0.txt")


plot(-rep1$V1/1000, rep1$V2*10000, ylim = c(0, 330000), xlim=c(-start_year, 0), type = 'l', xlab = "Years before present (x1000)", 
     ylab = "Effective population size", lwd=lwd, cex.lab=1.5, xaxt='n', col = 'chartreuse4', cex.axis=1.3)   

lines(-akikiki$V1/1000, akikiki$V2*10000, lwd=lwd, col = 'gray')
lines(-akekee$V1/1000, akekee$V2*10000, lwd=lwd, col = 'gold1')
lines(-poouli$V1/1000, poouli$V2*10000, lwd=lwd, col = 'coral4')
lines(-rep1$V1/1000, rep1$V2*10000, lwd=lwd, col = 'chartreuse4')

axis(side=1, cex.axis =1.3, labels = seq(start_year,0, by=-100), at = seq(-start_year,0, by=100))
legend(x = "topleft", legend = c("'Alala", "'Akikiki", "'Akeke'e", "Po'ouli"), col = c("chartreuse4", "gray", "gold1", "coral4"), lty = 1, bty = "n", lwd=5, cex = 1.5)

## plot replicates
for(rep in seq_along(alala_files)){
  rep_data <- alala_data[alala_data$rep == rep, ]
  #lines(-rep_data$V1/1000, rep_data$V2*10000, lwd = lwd, col = 'chartreuse4')
}

akikiki <- read.table("~/Documents/SDZWA/honeycreeper_genomics/analysis/PSMC/data/70AKIK_FilterB_Round2_26scaffolds_PASS_sortedbyRef_AKIK_1821-21127_bootstrap.0.txt")
akekee <- read.table("~/Documents/SDZWA/honeycreeper_genomics/analysis/PSMC/data/17AKEK_FilterB_Round3_29scaffolds_PASS_sortedbyRef_AKEK_1821-21146_bootstrap.0.txt")
poouli <- read.table("~/Documents/SDZWA/honeycreeper_genomics/analysis/PSMC/data/1Poouli_FilterB_Round1_30scaffolds_PASS_sortedbyRef_poouli_refsamp_bootstrap.0.txt")



dev.off()










last_gen <- alala_data[alala_data$V1==0,]
last_gen$V2*10000

rep1$V2*10


alala_last500k <- rep1[which(rep1$V1<500000 & rep1$V1>10000),]
mean(alala_last500k$V2)


