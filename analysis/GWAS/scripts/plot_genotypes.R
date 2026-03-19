

DLG1 <- c(0.56, 0.58,  0.0)

NEO1 <- c(0.603333333, 0.603333333, 0)


pdf("~/Documents/SDZWA/alala_genomics/analysis/GWAS/plots/recessive_lethal_genotypes.pdf",height=4, width=8 )
par(mfrow=c(1,2))
barplot(DLG1, ylab = "Proportion hatching", names.arg = c("hom ref", "het", "hom alt"), ylim=c(0,0.6), 
        main = "DLG1", col = 'blue4')

barplot(NEO1, ylab = "Proportion hatching", names.arg = c("hom ref", "het", "hom alt"), ylim=c(0,0.6), 
        main = "NEO1", col = 'blue4')
dev.off()