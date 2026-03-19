
setwd("~/Documents/SDZWA/alala_genomics/analysis/bird_heterozygosity/")
het_data <- read.csv("bird_heterozygosity.csv", header=T, stringsAsFactors=FALSE)
het_data <- het_data[order(het_data$pi, decreasing = T), ]




pdf("bird_species_heterozygosity.pdf", width = 6, height=11)
par(mar=c(5,14,2,2.5))

library(RColorBrewer)
cols <- brewer.pal(6, "YlOrRd")   
colors <- c("EX" = cols[6], "CR" = cols[5], "EN" = cols[4], "VU" = cols[3], "NT"= cols[2], "LC" = cols[1])
barplot(het_data$pi*1000, names.arg = het_data$Species, las=2, xlab = "Heterozygous sites per 1000 bp", col = colors[het_data$IUCN], horiz='T', 
        cex.axis = 1.3, cex.names = 1.3, cex.lab=1.5)
#legend("topright", legend = c("Extinct", "Critically endangered", "Endangered", "Vulnerable", "Near threatened", "Least concern" ), fill = colors, title = "IUCN status", cex=1.3, bty = "n")

dev.off()




pdf("bird_species_heterozygosity_presentation.pdf", width = 12, height=12)
par(mar=c(5,14,2,2.5))

library(RColorBrewer)
cols <- brewer.pal(6, "YlOrRd")   
colors <- c("EX" = cols[6], "CR" = cols[5], "EN" = cols[4], "VU" = cols[3], "NT"= cols[2], "LC" = cols[1], "EW" = 'black')
barplot(het_data$pi*1000, names.arg = het_data$Species, las=2, xlab = "Heterozygous sites per 1000 bp", col = colors[het_data$IUCN], horiz='T', 
        cex.axis = 1.3, cex.names = 1.3, cex.lab=1.5)

legend("topright", legend = c("Extinct", "Critically endangered", "Endangered", "Vulnerable", "Near threatened", "Least concern", "Extinct in the wild" ), 
       fill = colors, title = "IUCN status", cex=1.5, bty = "n")

dev.off()



