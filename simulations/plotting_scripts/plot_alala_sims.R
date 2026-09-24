
library(data.table)

sim_name <- "alala_noRelease_torgDFE_043026"
setwd(file.path("~/Documents/SDZWA/alala_genomics/simulations/data/", sim_name))

datafiles <- list.files(pattern = ".o1")

results_list <- vector("list", length(datafiles))

for (rep in seq_along(datafiles)) {
  
  # fast read
  data <- fread(datafiles[rep], skip = 20021)
  
  # find first index where popSizep2 > 0
  idx <- which(data$popSizep2 > 0)
  
  if (length(idx) == 0) next
  
  first <- idx[1]
  
  # keep only rows after first nonzero
  data_sub <- data[first:.N]
  
  # generate year vector
  data_sub[, year := 1973 + seq_len(.N)]
  
  # add rep column
  data_sub[, rep := rep]
  results_list[[rep]] <- data_sub
}

# combine once (fast)
data_frame <- rbindlist(results_list, fill = TRUE)

### subset data to match empirical trends
data_frame_keep <- data_frame[
  , if (any(year == 2005 & popSizep2 > 20 & popSizep2 < 100)) .SD,
  by = rep
]


# ensure it's a data.table
setDT(data_frame_keep)

# identify numeric columns only
num_cols <- names(data_frame_keep)[sapply(data_frame_keep, is.numeric)]

# optionally remove grouping variables from summary
num_cols <- setdiff(num_cols, c("year", "rep"))


means <- data_frame_keep[
  , lapply(.SD, mean, na.rm = TRUE),
  by = year,
  .SDcols = num_cols
]

sd <- data_frame_keep[
  , lapply(.SD, sd, na.rm = TRUE),
  by = year,
  .SDcols = num_cols
]

year <- 2100

#means$popSizep1[means$year>2030 & means$year<year] <- NA
means$FROH_1Mbp1[means$year>2005 & means$year<year] <- NA
means$B_genp1[means$year>2005 & means$year<year] <- NA
means$meanFitnessp1[means$year>2005 & means$year<year] <- NA
means$maxLetp1[means$year>2005 & means$year<year] <- NA

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
census_size_captive <- c(8, 10, 13, 18, 19, 21, 30, 38, 49, 51, 113, 120)
census_years_wild <- c(1978, 1983, 1992, 1993, 1994, 1996, 1998, 2000, 2002, 2004)
census_size_wild <- c(77, 25, 12, 12, 10, 5, 4, 4, 2, 0)
hatch_year <- c(1995,2005,2015)
FROH_1Mb <- c(0.286589952, 0.283749635, 0.326358927)

par(mfrow=c(4,1),mar = c(4,4.5,1,1), bty = "n")
### plot pop size
plot(means$year, means$popSizep2, type = "l", ylim = c(0,300), ylab = "Population size", lwd = mean_lwd, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year), col = col_mean_capt, xaxt='n')
#lines(means$year, means$deathsP1, col=col_mean_wild, lty=2)
#lines(means$year, means$deathsP2, col=col_mean_capt, lty=2)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$popSizep2 + sd$popSizep2, rev(means$popSizep2 - sd$popSizep2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)
#legend(x = "topleft", legend = c("Wild population", "Captive population"), col = c("blue4", "red4"), lty = 1, bty = "n", lwd=3, cex = 1)

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
# Pre-extinction wild population
pre_ext <- means$year <= 2005  # adjust as needed
polygon(c(means$year[pre_ext], rev(means$year[pre_ext])),
        c(means$FROH_1Mbp1[pre_ext] + sd$FROH_1Mbp1[pre_ext], 
          rev(means$FROH_1Mbp1[pre_ext] - sd$FROH_1Mbp1[pre_ext])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

# Post-release wild population  
post_rel <- means$year >= 2026  # adjust as needed
polygon(c(means$year[post_rel], rev(means$year[post_rel])),
        c(means$FROH_1Mbp1[post_rel] + sd$FROH_1Mbp1[post_rel], 
          rev(means$FROH_1Mbp1[post_rel] - sd$FROH_1Mbp1[post_rel])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)
points(hatch_year, FROH_1Mb, pch=19, col = col_mean_capt)

### plot fitness
plot(means$year, means$meanFitnessp2, type = "l", ylim = c(0.92,1.0), ylab = "Fitness", lwd = mean_lwd, cex.lab=1.2, xlab = "",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$meanFitnessp2 + sd$meanFitnessp2, rev(means$meanFitnessp2 - sd$meanFitnessp2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

lines(means$year, means$meanFitnessp2, lwd=mean_lwd, col = col_mean_capt)
lines(means$year, means$meanFitnessp1, lwd=mean_lwd, col = col_mean_wild)
# Pre-extinction wild population
pre_ext <- means$year <= 2005  # adjust as needed
polygon(c(means$year[pre_ext], rev(means$year[pre_ext])),
        c(means$meanFitnessp1[pre_ext] + sd$meanFitnessp1[pre_ext], 
          rev(means$meanFitnessp1[pre_ext] - sd$meanFitnessp1[pre_ext])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

# Post-release wild population  
post_rel <- means$year >= 2026  # adjust as needed
polygon(c(means$year[post_rel], rev(means$year[post_rel])),
        c(means$meanFitnessp1[post_rel] + sd$meanFitnessp1[post_rel], 
          rev(means$meanFitnessp1[post_rel] - sd$meanFitnessp1[post_rel])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

### plot inbreeding load 
plot(means$year, means$B_genp2, type = "l", ylim = c(0,4), ylab = "Inbreeding load", lwd = mean_lwd, cex.lab=1.2, xlab = "Year",xlim = c(start_year, end_year),  xaxt="n", col = col_mean_capt)
axis(1, at=seq(start_year,end_year, by=10), labels=c(seq(start_year, end_year, by = 10)))
polygon(c(means$year, rev(means$year)),
        c(means$B_genp2 + sd$B_genp2, rev(means$B_genp2 - sd$B_genp2)),
        col = adjustcolor(col_rep_capt, alpha.f = 0.5), border = NA)

lines(means$year, means$B_genp2, lwd=mean_lwd, col = col_mean_capt)
lines(means$year, means$B_genp1, lwd=mean_lwd, col = col_mean_wild)
# Pre-extinction wild population
pre_ext <- means$year <= 2005  # adjust as needed
polygon(c(means$year[pre_ext], rev(means$year[pre_ext])),
        c(means$B_genp1[pre_ext] + sd$B_genp1[pre_ext], 
          rev(means$B_genp1[pre_ext] - sd$B_genp1[pre_ext])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

# Post-release wild population  
post_rel <- means$year >= 2026  # adjust as needed
polygon(c(means$year[post_rel], rev(means$year[post_rel])),
        c(means$B_genp1[post_rel] + sd$B_genp1[post_rel], 
          rev(means$B_genp1[post_rel] - sd$B_genp1[post_rel])),
        col = adjustcolor(col_rep_wild, alpha.f = 0.5), border = NA)

dev.off()




means$popSizep1[means$year==2050]
median(data_frame_keep$popSizep1[data_frame_keep$year==2100])
hist(data_frame$popSizep1[data_frame$year==2100])


means$popSizep2[means$year==2036]
median(data_frame_keep$popSizep2[data_frame_keep$year==2036])







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











means$deathsP1[1:52]
means$birthsP1[1:52]
sum(means$deathsP1[means$year>2026 & means$year<2037])
sum(means$birthsP1[means$year>2026 & means$year<2037])

#par(mfrow=c(2,1))
#plot(means$year, means$birthsP2, xlim=c(1970, 2030), type="l", lwd=3, ylim=c(0,35))
#plot(means$year, means$deathsP2, xlim=c(1970, 2030), type="l", lwd=3, ylim=c(0,35))


pdf(paste("~/Documents/SDZWA/alala_genomics/simulations/plots/births_deaths_",sim_name,".pdf", sep=''), width = 5, height=5)
plot(means$year, means$popSizep1, xlim=c(2025, 2050), type="l", lwd=3, ylim=c(0,200), col = "black", xlab = "Year", ylab = "Population size", cex.lab=1.3)
lines(means$year, means$deathsP1, lwd=3, col = "red")
lines(means$year, means$birthsP1, lwd=3, col = "blue")
legend("topleft",
       legend = c("Population size", "Births", "Deaths"),
       col = c("black", "blue", "red"),
       bty = "n", lty=1, lwd=3, cex=1.2)

dev.off()




### plot for grant proposal

pdf("~/Documents/SDZWA/alala_genomics/simulations/plots/max_recessive_lethal_freq_presentation.pdf", height=5, width=4)
par(mfrow=c(2,1),mar = c(4,4.5,1,1), bty = "n")

end_year <- 2050

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

means$maxLetp2[means$year==2020]

dev.off()

