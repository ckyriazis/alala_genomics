library(dplyr)
library(ggplot2)

roh_combined_filtered <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/roh_combined_filtered.csv")
roh_combined_filtered_subset <- roh_combined_filtered %>% filter(grepl("K", individual_ID))

# pick chromosomes
chroms <- c("NC_063222.1", "NC_063225.1")

roh_plotdata <- roh_combined_filtered_subset %>%
  filter(chromosome %in% chroms,
         length > 1e5)   # only ROH >1 Mb

# get chromosome lengths
chr_lengths <- roh_plotdata %>%
  group_by(chromosome) %>%
  summarise(chr_len = max(ROH_end, na.rm = TRUE)) %>%
  arrange(chromosome)

# cumulative offsets
chr_lengths <- chr_lengths %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))

# join offsets
roh_plotdata <- roh_plotdata %>%
  left_join(chr_lengths, by = "chromosome") %>%
  mutate(
    ROH_start_global = ROH_start + chr_start,
    ROH_end_global   = ROH_end   + chr_start
  )

# ---- Heatmap prep ----
bin_size <- 5e5

bins <- chr_lengths %>%
  rowwise() %>%
  do({
    chr <- .$chromosome
    start <- .$chr_start
    len <- .$chr_len
    data.frame(
      chromosome = chr,
      bin_start = seq(start, start + len, by = bin_size),
      bin_end   = seq(start, start + len, by = bin_size) + bin_size
    )
  }) %>%
  ungroup()

roh_bins <- bins %>%
  rowwise() %>%
  mutate(
    count = sum(roh_plotdata$ROH_start_global <= bin_end &
                  roh_plotdata$ROH_end_global   >= bin_start, na.rm = TRUE)
  )

# sort individuals by numeric ID
roh_plotdata$ID_num <- as.numeric(gsub("[^0-9]", "", roh_plotdata$individual_ID))
roh_plotdata <- roh_plotdata[order(roh_plotdata$ID_num), ]
roh_plotdata$individual_ID <- factor(
  roh_plotdata$individual_ID,
  levels = rev(unique(roh_plotdata$individual_ID))
)

# ---- Gene marker prep ----
# NC_063222.1 = chr10 (DLG1), NC_063225.1 = chr13 (NEO1)
gene_markers <- data.frame(
  gene       = c("DLG1", "NEO1"),
  chromosome = c("NC_063222.1", "NC_063225.1"),
  position   = c(15792527, 20914015)
) %>%
  left_join(chr_lengths, by = "chromosome") %>%
  mutate(global_pos = position + chr_start)

# ---- Genotype data (from metadata CSV) ----
# 0 = hom ref, 1 = het, 2 = hom alt (recessive lethal), NA = missing
metadata <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/alala_medata.csv")
genotypes <- metadata %>%
  select(individual_ID = Sample, DLG1, NEO1)

# Build long-format genotype data for points on the plot
# Only include individuals present in roh_plotdata
geno_long <- genotypes %>%
  filter(individual_ID %in% levels(roh_plotdata$individual_ID)) %>%
  tidyr::pivot_longer(cols = c(DLG1, NEO1), names_to = "gene", values_to = "genotype") %>%
  left_join(gene_markers %>% select(gene, global_pos), by = "gene") %>%
  left_join(
    data.frame(
      individual_ID = levels(roh_plotdata$individual_ID),
      y_pos = seq_along(levels(roh_plotdata$individual_ID))
    ),
    by = "individual_ID"
  ) %>%
  filter(!is.na(genotype), !is.na(global_pos)) %>%
  mutate(genotype_label = case_when(
    genotype == 0 ~ "hom ref (0)",
    genotype == 1 ~ "het (1)",
    genotype == 2 ~ "hom alt (2)"
  ))

# ---- Plot ----
setwd("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/")
pdf("alala_ROH_plot_embryos_with_genotypes.pdf", width = 8, height = 9)

ggplot() +
  # heatmap track
  geom_tile(data = roh_bins,
            aes(x = (bin_start + bin_end) / 2,
                y = length(unique(roh_combined_filtered_subset$individual_ID)) + 3,
                fill = count,
                width = bin_size, height = 2)) +
  scale_fill_viridis_c(name = "ROH count") +
  
  # individual ROH tracks
  geom_segment(data = roh_plotdata,
               aes(x = ROH_start_global, xend = ROH_end_global,
                   y = as.numeric(individual_ID),
                   yend = as.numeric(individual_ID)),
               color = "steelblue", size = 2) +
  
  # vertical chromosome boundaries
  geom_vline(xintercept = chr_lengths$chr_start, color = "black", lwd = 0.5) +
  
  # genotype points per individual at gene positions (circles for both genes)
  ggnewscale::new_scale_fill() +
  geom_point(data = geno_long,
             aes(x = global_pos, y = y_pos,
                 fill = genotype_label),
             shape = 21, size = 2, stroke = 0.3, color = "grey30", alpha = 0.6) +
  scale_fill_manual(
    name = "Genotype",
    values = c("hom ref (0)" = "white",
               "het (1)"     = "orange",
               "hom alt (2)" = "red"),
    na.value = "transparent",
    guide = guide_legend(override.aes = list(shape = 21, size = 3, color = "grey30"))
  ) +
  
  scale_y_continuous(breaks = 1:length(levels(roh_plotdata$individual_ID)),
                     labels = levels(roh_plotdata$individual_ID)) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(x = "Chromosome",
       y = "Individual",
       title = "")

dev.off()