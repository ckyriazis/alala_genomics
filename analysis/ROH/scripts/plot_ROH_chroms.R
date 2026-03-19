
# Load required package
library(dplyr)
library(ggplot2)

# Set the directory containing your ROH files
roh_dir <- "~/Documents/SDZWA/alala_genomics/analysis/ROH/data/175AL_round1/"

# List all files (adjust pattern if needed, e.g. ".txt" or ".csv")
roh_files <- list.files(roh_dir, full.names = TRUE, pattern = "AL")

read_roh <- function(file) {
  # read with fill=TRUE so ragged rows get padded
  df <- read.table(gzfile(file), header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
  
  # make sure df has at least 8 columns
  if (ncol(df) < 8) {
    df <- cbind(df, matrix(NA, nrow(df), 8 - ncol(df)))
  }
  
  # assign 8 column names
  colnames(df)[1:8] <- c("row_type", "sample", "chrom", "start", "end",
                         "length", "num_markers", "qual")
  
  df %>%
    filter(row_type == "RG") %>%
    select(individual_ID = sample,
           chromosome   = chrom,
           ROH_start    = start,
           ROH_end      = end,
           length = length)
}

roh_combined <- roh_files %>%
  lapply(read_roh) %>%
  bind_rows()

roh_combined_filtered <- roh_combined %>% filter(length>1e6)

#write.csv(roh_combined_filtered, file="~/Documents/SDZWA/alala_genomics/analysis/ROH/roh_combined_filtered.csv")

roh_combined_filtered <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/roh_combined_filtered_merged.csv")
head(roh_combined_filtered)

roh_combined_filtered_subset <- roh_combined_filtered %>% filter(grepl("AL", individual_ID))

hist(log(roh_combined_filtered$length), breaks =100)


# pick chromosomes
chroms <- unique(roh_combined_filtered_subset$chromosome)[1:5]

roh_plotdata <- roh_combined_filtered_subset %>%
  filter(chromosome %in% chroms,
         length > 1e6)   # only ROH >1 Mb

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
# bin size (adjustable)
bin_size <- 5e5

# make bins for each chromosome
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

# count how many individuals overlap each bin
roh_bins <- bins %>%
  rowwise() %>%
  mutate(
    count = sum(roh_plotdata$ROH_start_global <= bin_end &
                  roh_plotdata$ROH_end_global   >= bin_start, na.rm = TRUE)
  )

# extract numeric part of ID using gsub
roh_plotdata$ID_num <- as.numeric(gsub("[^0-9]", "", roh_plotdata$individual_ID))

# sort by numeric ID
roh_plotdata <- roh_plotdata[order(roh_plotdata$ID_num), ]

# reset factor levels so ggplot respects the order
roh_plotdata$individual_ID <- factor(
  roh_plotdata$individual_ID,
  levels = unique(roh_plotdata$individual_ID)
)


setwd("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/")
pdf("alala_ROH_plot_adults.pdf", width=8, height=9)

# ---- Plot ----
roh_plotdata$individual_ID <- factor(
  roh_plotdata$individual_ID,
  levels = rev(unique(roh_plotdata$individual_ID))
)

ggplot() +
  # heatmap track (population-level ROH abundance)
  geom_tile(data = roh_bins, aes(x = (bin_start + bin_end) / 2,
                                 y = length(unique(roh_combined_filtered_subset$individual_ID))+3,  # place heatmap at top
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
  geom_vline(xintercept = chr_lengths$chr_start, color = "black", lwd=0.5) +
  
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







  