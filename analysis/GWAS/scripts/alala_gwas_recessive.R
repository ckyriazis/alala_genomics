# Recessive Model GWAS for Embryo Mortality in 'Alala
# Testing for homozygous derived alleles causing mortality

# Required packages
library(SNPRelate)
library(qqman)
library(data.table)
library(gdsfmt)
library(ggplot2)
library(dplyr)

setwd("~/Documents/SDZWA/alala_genomics/analysis/GWAS/data/")

# ============================================================================
# Load data
# ============================================================================

plink_prefix <- "175AL_autosomes_PASS_SNPs_noFilter"

# --- Load BIM (always as character, avoid factor issues) ---
bim <- read.table(
  paste0(plink_prefix, ".bim"),
  header = FALSE, stringsAsFactors = FALSE,
  col.names = c("chr", "snp", "cm", "pos", "A1", "A2")
)

# --- Save original BIM without modification ---
orig_bim <- bim

# --- Build a stable mapping: sort chromosome names consistently ---
chr_unique <- sort(unique(bim$chr))

chr_map <- data.frame(
  chr = chr_unique,
  chr.num  = seq_along(chr_unique)
)

cat("\nChromosome mapping:\n")
print(chr_map)

# --- Create temporary BIM for GDS conversion ---
bim_tmp <- bim

# Convert chromosome names → integers
bim_tmp$chr <- match(bim_tmp$chr, chr_unique)

# Assign stable integer SNP IDs
bim_tmp$snp <- seq_len(nrow(bim_tmp))

tmp_bim_path <- paste0(plink_prefix, ".bim.tmp")

write.table(
  bim_tmp, tmp_bim_path,
  row.names = FALSE, col.names = FALSE,
  quote = FALSE, sep = "\t"
)

# --- Convert to GDS using the temp BIM ---
snpgdsBED2GDS(
  bed.fn = paste0(plink_prefix, ".bed"),
  bim.fn = tmp_bim_path,
  fam.fn = paste0(plink_prefix, ".fam"),
  out.gdsfn = "alala.gds"
)

# --- Remove temporary BIM (original BIM stays untouched) ---
file.remove(tmp_bim_path)

# --- Open GDS ---
genofile <- snpgdsOpen("alala.gds")

# --- Load phenotype data ---
pheno <- read.csv("samples_meta.csv", stringsAsFactors = FALSE)

sample_ids_gds <- read.gdsn(index.gdsn(genofile, "sample.id"))
pheno_ordered <- pheno[match(sample_ids_gds, pheno$sample.id), ]

# --- Create consistent SNP info table using original BIM ---
snp_info <- orig_bim
snp_info$snpID <- seq_len(nrow(snp_info))

cat("\n=== Recessive Model GWAS ===\n")
cat("Total samples:", nrow(pheno_ordered), "\n")
cat("Cases (died in egg):", sum(pheno_ordered$embryo_mortality == 0, na.rm = TRUE), "\n")
cat("Controls (hatched):", sum(pheno_ordered$embryo_mortality == 1, na.rm = TRUE), "\n")
cat("Total SNPs:", nrow(snp_info), "\n\n")

# ============================================================================
# STEP 1: Extract all genotypes at once (more efficient)
# ============================================================================

cat("Extracting genotypes...\n")
all_genos <- snpgdsGetGeno(genofile, with.id = TRUE)

# Verify sample order matches
if(!identical(all_genos$sample.id, sample_ids_gds)) {
  stop("Sample order mismatch!")
}

# ============================================================================
# STEP 2: Recessive model testing for each SNP
# ============================================================================

cat("Running recessive model tests...\n")

# Initialize results
n_snps <- nrow(snp_info)
results <- data.frame(
  snpID = snp_info$snpID,
  chr = snp_info$chr,
  pos = snp_info$pos,
  A1 = snp_info$A1,
  A2 = snp_info$A2,
  n_total = NA,
  n_died = NA,
  n_hatched = NA,
  hom_alt_died = NA,
  hom_alt_hatched = NA,
  het_died = NA,
  het_hatched = NA,
  hom_ref_died = NA,
  hom_ref_hatched = NA,
  freq_hom_alt_died = NA,
  freq_hom_alt_hatched = NA,
  fisher_p = NA,
  OR = NA,
  stringsAsFactors = FALSE
)

# Progress bar
pb <- txtProgressBar(min = 0, max = n_snps, style = 3)

for(i in 1:n_snps) {
  # Get genotypes for this SNP (all samples for SNP i)
  # Genotype matrix is [samples x SNPs], so we need column i
  genos <- all_genos$genotype[, i]
  
  # Create data frame
  geno_df <- data.frame(
    genotype = genos,
    phenotype = pheno_ordered$embryo_mortality
  )
  
  # Remove missing data
  geno_df <- geno_df[!is.na(geno_df$genotype) & !is.na(geno_df$phenotype), ]
  
  if(nrow(geno_df) < 100) {
    next  # Skip SNPs with too much missing data
  }
  
  # Count genotypes by phenotype
  # 0 = died in egg, 1 = hatched
  # genotype: 0 = hom ref, 1 = het, 2 = hom alt
  
  died <- geno_df$phenotype == 0
  hatched <- geno_df$phenotype == 1
  
  hom_alt_died <- sum(geno_df$genotype == 2 & died)
  hom_alt_hatched <- sum(geno_df$genotype == 2 & hatched)
  het_died <- sum(geno_df$genotype == 1 & died)
  het_hatched <- sum(geno_df$genotype == 1 & hatched)
  hom_ref_died <- sum(geno_df$genotype == 0 & died)
  hom_ref_hatched <- sum(geno_df$genotype == 0 & hatched)
  
  n_died_total <- sum(died)
  n_hatched_total <- sum(hatched)
  
  # FILTER: Only test if hom alt is MORE frequent in died than hatched
  # This is the key filter for recessive lethal direction
  freq_hom_alt_died_temp <- hom_alt_died / n_died_total
  freq_hom_alt_hatched_temp <- hom_alt_hatched / n_hatched_total
  
  if(freq_hom_alt_hatched_temp > 0.05 | freq_hom_alt_died_temp < 0.05) {
    next  # Skip SNPs with more than x hom alt individuals that hatched
  }
  
  # Store counts
  results$n_total[i] <- nrow(geno_df)
  results$n_died[i] <- n_died_total
  results$n_hatched[i] <- n_hatched_total
  results$hom_alt_died[i] <- hom_alt_died
  results$hom_alt_hatched[i] <- hom_alt_hatched
  results$het_died[i] <- het_died
  results$het_hatched[i] <- het_hatched
  results$hom_ref_died[i] <- hom_ref_died
  results$hom_ref_hatched[i] <- hom_ref_hatched
  
  # Calculate frequencies
  results$freq_hom_alt_died[i] <- freq_hom_alt_died_temp
  results$freq_hom_alt_hatched[i] <- freq_hom_alt_hatched_temp
  
  # Recessive model: compare hom alt vs (het + hom ref)
  # Create 2x2 contingency table
  contingency <- matrix(c(
    hom_alt_died, het_died + hom_ref_died,      # died
    hom_alt_hatched, het_hatched + hom_ref_hatched  # hatched
  ), nrow = 2, byrow = TRUE)
  
  # Fisher's exact test (better for small sample sizes)
  if(all(contingency >= 0) && sum(contingency) > 0) {
    fisher_result <- fisher.test(contingency)
    results$fisher_p[i] <- fisher_result$p.value
    results$OR[i] <- fisher_result$estimate
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)
# ============================================================================
# STEP 3: Filter and rank results
# ============================================================================

cat("\n\nFiltering results...\n")

# Remove SNPs with missing p-values
results_clean <- results[!is.na(results$fisher_p), ]
 
# Sort by p-value
results_sorted <- results_clean[order(results_clean$fisher_p), ]
results_sorted <- read.table("~/Documents/SDZWA/alala_genomics/analysis/GWAS/output/recessive_gwas_results_full_nofilter0.05.txt", header = T)

# Add -log10 p-value
results_sorted$log10p <- -log10(results_sorted$fisher_p)

# Bonferroni correction
bonf_threshold <- 0.05 / nrow(results_sorted)
suggestive_threshold <- 1e-5  # More lenient for this specific test

# Significant at 5% FDR
results_sorted$fdr <- p.adjust(results_sorted$fisher_p, method = "BH")
fdr_threshold <- 0.05
fdr_sig <- results_sorted[results_sorted$fdr < fdr_threshold, ]


cat("\nTotal SNPs tested:", nrow(results_sorted), "\n")
cat("Bonferroni threshold:", format(bonf_threshold, scientific = TRUE), "\n")
cat("Suggestive threshold:", format(suggestive_threshold, scientific = TRUE), "\n")
cat("Bonferroni significant SNPs:", sum(results_sorted$fisher_p < bonf_threshold), "\n")
cat("FDR 5% significant:", sum(results_sorted$fdr < fdr_threshold), "\n")
cat("Suggestive SNPs:", sum(results_sorted$fisher_p < suggestive_threshold), "\n\n")


# ============================================================================
# STEP 4: Save results
# ============================================================================

# Save all results
write.table(results_sorted, 
            "~/Documents/SDZWA/alala_genomics/analysis/GWAS/output/recessive_gwas_results_full_nofilter0.05.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Save top 500 for visualization
top_500 <- head(results_sorted, 500)
write.table(top_500,
            "~/Documents/SDZWA/alala_genomics/analysis/GWAS/output/recessive_gwas_top500_nofilter0.05.txt",
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# ============================================================================
# STEP 5: Manhattan plot
# ============================================================================

# Create chromosome mapping for plotting
chr_unique <- sort(unique(results_sorted$chr))
chr_map <- data.frame(
  chr = chr_unique,
  chr.num = 1:length(chr_unique)
)

manhattan_data <- merge(results_sorted, chr_map, by.x = "chr", by.y = "chr")
manhattan_data <- data.frame(
  SNP = manhattan_data$snpID,
  CHR = manhattan_data$chr.num,
  BP = manhattan_data$pos,
  P = manhattan_data$fisher_p
)

manhattan_data <- manhattan_data[!is.na(manhattan_data$P) & manhattan_data$BP > 0, ]

# Filter to only first x chromosomes
manhattan_data_filtered <- manhattan_data %>%
  filter(CHR <= 30)

# Prepare data for ggplot Manhattan plot
manhattan_data_plot <- manhattan_data_filtered %>%
  # Calculate cumulative position for x-axis
  group_by(CHR) %>%
  summarise(chr_len = max(BP)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(manhattan_data_filtered, ., by = "CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot,
         log10p = -log10(P))

# Calculate axis positions (center of each chromosome)
# Modify axis_df to show only odd chromosomes
axis_df <- manhattan_data_plot %>%
  group_by(CHR) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2) %>%
  mutate(label = ifelse(CHR %% 2 == 0 | CHR<16, as.character(CHR), ""))

# Create the plot
p <- ggplot(manhattan_data_plot, aes(x = BPcum, y = log10p)) +
  # Points colored by chromosome
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.5) +
  scale_color_manual(values = rep(c("gray60", "lightblue3"), 15)) +  # 26/2 = 13 pairs
  
  # Threshold lines
  geom_hline(yintercept = -log10(suggestive_threshold), 
             linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_hline(yintercept = -log10(bonf_threshold), 
             linetype = "dashed", color = "red", linewidth = 0.5) +
  
  # Customize axes
  scale_x_continuous(label = axis_df$label, breaks = axis_df$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,6.5)) +
  
  # Labels and theme
  labs(x = "Chromosome", 
       y = expression(-log[10](italic(p)))) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

# Display the plot
print(p)

# Save the plot
ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/plots/recessive_manhattan_plot_nofilter_freq0.05.png", 
       plot = p, 
       width = 10, 
       height = 3.5, 
       dpi = 300)


# ============================================================================
# STEP 6: Visualize chrom 10 region & LD
# ============================================================================

# Get all chr10 SNPs with any signal
chr10_results <- results_sorted %>%
  filter(chr == "NC_063222.1",
         pos >= 14000000,
         pos <= 24000000) %>%
  arrange(pos)

# Plot p-value across region
ggplot(chr10_results, aes(x = pos/1e6, y = -log10(fisher_p))) +
  geom_point(aes(color = -log10(fisher_p) > 4), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("gray70", "darkred"), guide = "none") +
  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue") +
  labs(title = "",
       x = "Position (Mb)",
       y = "-log10(p-value)") +
  theme_minimal(base_size = 14) +
  annotate("rect", xmin = 15.7, xmax = 15.85, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "red") +
  #annotate("rect", xmin = 22.4, xmax = 22.5, ymin = 0, ymax = Inf, 
  #         alpha = 0.1, fill = "red") +
  annotate("text", x = 15.8, y = 5.5, label = "DLG1", size = 3) 
  #annotate("text", x = 22.45, y = 5.5, label = "Region 2", size = 3)

ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/plots/chr10_zoomed.png", width = 10, height = 4, dpi = 300)



# Pick reference SNP (DLG1 missense)
dlg1_snp_id <- which(snp_info$chr == "NC_063222.1" & 
                       snp_info$pos == 15792527)
dlg1_genos <- all_genos$genotype[, dlg1_snp_id]

# Get all SNPs in chr10 region (sample every 50th for speed)
chr10_snp_ids <- snp_info %>%
  filter(chr == "NC_063222.1",
         pos >= 14000000,
         pos <= 24000000) %>%
  slice(seq(1, n(), by = 5))

# Calculate LD with DLG1 for each SNP
ld_decay <- sapply(chr10_snp_ids$snpID, function(snp_id) {
  genos <- all_genos$genotype[, snp_id]
  r <- cor(dlg1_genos, genos, use = "complete.obs")
  r^2
})

# Create LD decay plot
ld_decay_df <- data.frame(
  pos = chr10_snp_ids$pos / 1e6,
  r2 = ld_decay,
  distance_from_dlg1 = abs(chr10_snp_ids$pos - 15792527) / 1e6
)

# Plot LD decay
ggplot(ld_decay_df, aes(x = distance_from_dlg1, y = r2)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "loess", span = 0.3, se = TRUE, color = "red") +
  labs(title = "LD decay from DLG1 p.Cys23Tyr",
       x = "Distance from DLG1 (Mb)",
       y = "r² with DLG1") +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1))

# Also plot by pos
ggplot(ld_decay_df, aes(x = pos, y = r2)) +
  geom_point(aes(color = r2), alpha = 0.7, size = 2) +
  geom_vline(xintercept = 15.79, linetype = "dashed", color = "red") +
  scale_color_gradient2(low = "gray80", mid = "orange", high = "darkred",
                        midpoint = 0.5) +
  labs(title = "",
       x = "Position (Mb)",
       y = "r² with DLG1") +
  theme_minimal(base_size = 14)

ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/plots/chr10_ld_with_dlg1.png", width = 10, height = 4, dpi = 300)


### plot heterozygosity as a test of inversion 

# Get individuals by genotype at DLG1
dlg1_snp_id <- which(snp_info$chr == "10" & 
                       snp_info$pos == 15792527)
dlg1_genos <- all_genos$genotype[, dlg1_snp_id]

# Classify individuals
hom_ref <- which(dlg1_genos == 0)      # Standard/Standard
het <- which(dlg1_genos == 1)           # Standard/Inverted
hom_alt <- which(dlg1_genos == 2)       # Inverted/Inverted

# Get all SNPs across the region
region_snps <- snp_info %>%
  filter(chr == "10",
         pos >= 14000000,
         pos <= 24000000)

# Calculate heterozygosity in each group
calc_het <- function(geno_ids) {
  region_genos <- all_genos$genotype[geno_ids, region_snps$snpID]
  het_rate <- apply(region_genos, 2, function(x) {
    sum(x == 1, na.rm = TRUE) / sum(!is.na(x))
  })
  data.frame(
    pos = region_snps$pos,
    het_rate = het_rate
  )
}

het_hom_ref <- calc_het(hom_ref) %>% mutate(group = "Hom Ref (Std/Std)")
het_het <- calc_het(het) %>% mutate(group = "Het (Std/Inv)")
het_hom_alt <- calc_het(hom_alt) %>% mutate(group = "Hom Alt (Inv/Inv)")

# Combine
het_all <- bind_rows(het_hom_ref, het_hom_alt,het_het)

# Plot
ggplot(het_all, aes(x = pos/1e6, y = het_rate, color = group)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_vline(xintercept = c(14.9, 15.8, 22.4, 22.5), 
             linetype = "dashed", alpha = 0.5) +
  labs(title = "Heterozygosity patterns across Chr10 region",
       subtitle = "Inversion signature: reduced Het in heterozygotes",
       x = "Position (Mb)",
       y = "Heterozygosity rate",
       color = "Genotype") +
  theme_minimal(base_size = 14)

ggsave("chr10_heterozygosity_by_genotype.png", width = 12, height = 6, dpi = 300)





# Get all SNPs across the region
region_snps <- snp_info %>%
  dplyr::filter(
    chr == "10",
    pos >= 8000000,
    pos <= 10000000
  )

pca_region <- snpgdsPCA(
  genofile,
  snp.id = region_snps$snpID,
  num.thread = 4
)

# Data frame
pca_df <- data.frame(
  sample_id = pca_region$sample.id,
  PC1 = pca_region$eigenvect[, 1],
  PC2 = pca_region$eigenvect[, 2],
  PC3 = pca_region$eigenvect[, 3],
  dlg1_genotype = dlg1_genos
)


ggplot(pca_df, aes(x = PC1, y = PC2, color = factor(dlg1_genotype))) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("blue", "purple", "red"),
                     labels = c("Std/Std", "Std/Inv", "Inv/Inv")) +
  labs(title = "PCA of Chr10 region",
       color = "Genotype") +
  theme_minimal(base_size = 14)







### LD heatmap 

region_snps <- snp_info %>%
  dplyr::filter(
    chr == "10",
    pos >= 10000000,
    pos <= 25000000
  )

snp_ids <- region_snps$snpID

window_snps <- region_snps %>%
  filter(pos >= 10000000 & pos <= 25000000)

ld <- snpgdsLDMat(genofile, snp.id = window_snps$snpID)

image(ld$LD^2)
ld_matrix <- ld$LD^2


ld_df <- reshape2::melt(ld_matrix)
colnames(ld_df) <- c("SNP1", "SNP2", "r2")

# Add genomic positions for nicer axis labels
ld_df$pos1 <- region_snps$pos[ld_df$SNP1]
ld_df$pos2 <- region_snps$pos[ld_df$SNP2]


ggplot(ld_df, aes(x = pos1, y = pos2, fill = r2)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = expression(r^2)) +
  labs(
    title = "LD Heatmap (Chr10 region)",
    x = "Genomic position",
    y = "Genomic position"
  ) +
  theme_minimal(base_size = 14) +
  coord_fixed()



