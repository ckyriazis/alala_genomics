
sim_name <- "alala_release10_0.1_10yr_020626"

setwd(paste("~/Documents/SDZWA/alala_genomics/simulations/data/", sim_name, sep=""))
datafiles <- list.files(pattern=".o1")

data_frame <- data.frame(matrix(nrow = 0, ncol = 21)) #15

for(rep in seq_along(datafiles)){
  data <- read.table(datafiles[rep], sep=",", header=T, skip=20021) #20020
  count <- 0
  year <- 1973
  for(row in 1:length(data$popSizep2)){
    if(data$popSizep2[row] > 0){
      count <- 1
      year <- year+1
    }
    if(count == 1){
      data_frame <- rbind(data_frame,cbind(data[row,],rep, year))
    }
  }
}

### subset data to match empirical trends
data_frame_keep <- data.frame(matrix(nrow = 0, ncol = 22)) #16

for(rep in seq_along(datafiles)){
  data <- data_frame[data_frame$rep == rep,]
  if(data$popSizep2[data$year==2005]>20 && data$popSizep2[data$year==2005]<100){
    data_frame_keep <- rbind(data_frame_keep, data)
  }
}


### get means
years <- unique(data_frame_keep$year)
means <- data.frame(matrix(nrow = 0, ncol = 22))
sd <- data.frame(matrix(nrow = 0, ncol = 22))

for(year in years){
  data_this_gen <- data_frame_keep[which(data_frame_keep$year==year),]
  means <- rbind(means,colMeans(data_this_gen, na.rm = T))
  sd <- rbind(sd, sapply(data_this_gen, sd, na.rm = T))
  
}
colnames(means) <- colnames(data_frame)
colnames(sd) <- colnames(data_frame)
sd$year <- means$year
sd$gen <- means$gen


year <- 2100

means$FROH_1Mbp1[means$year>1999 & means$year<year] <- NA
means$B_genp1[means$year>1999 & means$year<year] <- NA
means$meanFitnessp1[means$year>1999 & means$year<year] <- NA
means$maxLetp1[means$year>1999 & means$year<year] <- NA
#sd$FROH_1Mbp1[sd$year>2005 & sd$year<2025] <- NA
#sd$B_genp1[sd$year>2005 & sd$year<2025]<- NA
#sd$meanFitnessp1[sd$year>2005 & sd$year<2025]<- NA

### how many reps remain
length(unique(data_frame_keep$rep))
length(unique(data_frame$rep))



pdf(paste("~/Documents/SDZWA/alala_genomics/simulations/plots/plot_",sim_name,".pdf", sep=''), height=6, width=3)

start_year <- 1970
end_year <- 2100
col_mean_wild <- "blue4"
col_mean_capt <- "red4"
col_rep_wild <- "blue4"
col_rep_capt <- "red4"

sd <- sd[sd$year<=end_year,]
means <- means[means$year<=end_year,]


mean_lwd <- 2

census_years_captive <- c(1978, 1983, 1992, 1994, 1996, 1998, 2000, 2002, 2004, 2006,2014, 2025)
census_size_captive <- c(8, 10, 13, 18, 19, 21, 30, 38, 49, 51,113, 120)
census_years_wild <- c(1978, 1983, 1992, 1993, 1994, 1996, 1998, 2000, 2002, 2004)
census_size_wild <- c(77, 25, 12, 12, 10, 5, 4, 4, 2, 0)
hatch_year <- c(1995,2005,2015)
FROH_1Mb <- c(0.286589952, 0.283749635, 0.326358927)

par(mfrow=c(4,1),mar = c(4,4.5,1,1), bty = "n")
### plot pop size
plot(means$year, means$popSizep2, type = "l", ylim = c(0,300), ylab = "Population size", lwd = mean_lwd, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year), col = col_mean_capt, xaxt='n')
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$popSizep2 + sd$popSizep2, rev(means$popSizep2 - sd$popSizep2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)
legend(x = "topleft", legend = c("Wild population", "Breeding population"), col = c("blue4", "red4"), 
       lty = 1, bty = "n", lwd=3, cex = 1)

lines(means$year, means$popSizep1, lwd=mean_lwd, col = col_mean_wild)
polygon(c(means$year, rev(means$year)),
        c(means$popSizep1 + sd$popSizep1, rev(means$popSizep1 - sd$popSizep1)),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

points(census_years_wild, census_size_wild, pch=19, col = col_mean_wild)

lines(means$year, means$popSizep2, lwd=mean_lwd, col = col_mean_capt)
points(census_years_captive, census_size_captive, pch=19, col = col_mean_capt)

### plot FROH
plot(means$year, means$FROH_1Mbp2, type = "l", ylim = c(0,0.45), ylab = expression(' F'[ROH]), lwd = mean_lwd, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$FROH_1Mbp2 + sd$FROH_1Mbp2, rev(means$FROH_1Mbp2 - sd$FROH_1Mbp2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

lines(means$year, means$FROH_1Mbp1, lwd=mean_lwd, col = col_mean_wild)
#polygon(c(means$year, rev(means$year)),
#        c(means$FROH_1Mbp1 + sd$FROH_1Mbp1, rev(means$FROH_1Mbp1 - sd$FROH_1Mbp1)),
#        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)
points(hatch_year, FROH_1Mb, pch=19, col = col_mean_capt)

### plot fitness
plot(means$year, means$meanFitnessp2, type = "l", ylim = c(0.92,1.0), ylab = "Fitness", lwd = mean_lwd, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$meanFitnessp2 + sd$meanFitnessp2, rev(means$meanFitnessp2 - sd$meanFitnessp2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

lines(means$year, means$meanFitnessp2, lwd=mean_lwd, col = col_mean_capt)
lines(means$year, means$meanFitnessp1, lwd=mean_lwd, col = col_mean_wild)
#polygon(c(means$year, rev(means$year)),
#        c(means$meanFitnessp1 + sd$meanFitnessp1, rev(means$meanFitnessp1 - sd$meanFitnessp1)),
#        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

### plot inbreeding load 
plot(means$year, means$B_genp2, type = "l", ylim = c(0,4), ylab = "Inbreeding load", lwd = mean_lwd, cex.lab=1.2, xlab = "Year",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$B_genp2 + sd$B_genp2, rev(means$B_genp2 - sd$B_genp2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

lines(means$year, means$B_genp2, lwd=mean_lwd, col = col_mean_capt)
lines(means$year, means$B_genp1, lwd=mean_lwd, col = col_mean_wild)
#polygon(c(means$year, rev(means$year)),
#        c(means$B_genp1 + sd$B_genp1, rev(means$B_genp1 - sd$B_genp1)),
#        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)


dev.off()




means$popSizep1[means$year==2100]


means$FROH_1Mbp2[means$year==1975]
means$FROH_1Mbp2[means$year==2000]
means$FROH_1Mbp2[means$year==2025]
means$FROH_1Mbp2[means$year==2100]


means$meanFitnessp2[means$year==1975]
means$meanFitnessp2[means$year==2000]
means$meanFitnessp2[means$year==2025]
means$meanFitnessp2[means$year==2100]


means$B_genp2[means$year==1975]
means$B_genp2[means$year==2025]
means$B_genp2[means$year==2100]


means$maxLetp2[means$year==1975]
means$maxLetp2[means$year==2025]
means$maxLetp2[means$year==2100]


means$maxLetp2[means$year==1990]
means$maxLetp2[means$year==2020]



pdf("~/Documents/SDZWA/alala_genomics/simulations/plots/max_recessive_lethal_freq.pdf", height=2.5, width=6)
par(mar=c(4,4,1,1))

### plot maximum lethal allele frequency
plot(means$year, means$maxLetp2, type = "l", ylim = c(0,0.6), ylab = "Max lethal freq", lwd = 5, cex.lab=1.2, xlab = "Year",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$maxLetp2 + sd$maxLetp2*2, rev(means$maxLetp2 - sd$maxLetp2*2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

#for(file in seq_along(datafiles)){
#  data <- data_frame[which(data_frame$rep==file),]
#  lines(data$year, data$maxLetp2, col = col_rep_capt, lwd = 0.6)
#}

#lines(means$year, means$maxLetp2, lwd=10, col = col_rep_capt)
#lines(means$year, means$maxLetp1, lwd=4, col = col_mean_wild)
#polygon(c(means$year, rev(means$year)),
#        c(means$B_genp1 + sd$B_genp1, rev(means$B_genp1 - sd$B_genp1)),
#        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

dev.off()







pdf("~/Documents/SDZWA/alala_genomics/simulations/plots/max_recessive_lethal_freq_presentation.pdf", height=5, width=4)
par(mfrow=c(2,1),mar = c(4,4.5,1,1), bty = "n")


### plot pop size
plot(means$year, means$popSizep2, type = "l", ylim = c(0,130), ylab = "Population size", lwd = 5, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year), col = col_mean_capt, xaxt='n')
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$popSizep2 + sd$popSizep2, rev(means$popSizep2 - sd$popSizep2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)


### plot maximum lethal allele frequency
plot(means$year, means$maxLetp2, type = "l", ylim = c(0,0.6), ylab = "Max lethal freq", lwd = 5, cex.lab=1.2, xlab = "Year",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$maxLetp2 + sd$maxLetp2*2, rev(means$maxLetp2 - sd$maxLetp2*2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)


dev.off()




