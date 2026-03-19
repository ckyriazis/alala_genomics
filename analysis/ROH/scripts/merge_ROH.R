library(dplyr)
library(tidyr)

# Read in your ROH table
roh <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/roh_combined_filtered.csv")

# Ensure data is ordered by individual, chromosome, and start position
roh <- roh %>%
  arrange(individual_ID, chromosome, ROH_start)

# Define maximum gap (in base pairs)
max_gap <- 500000

# Merge ROHs that are close together
merged_roh <- roh %>%
  group_by(individual_ID, chromosome) %>%
  # Assign group IDs where gaps > max_gap cause a new group
  mutate(group_id = cumsum(c(TRUE, (ROH_start[-1] - ROH_end[-n()]) > max_gap))) %>%
  group_by(individual_ID, chromosome, group_id) %>%
  summarise(
    ROH_start = min(ROH_start),
    ROH_end   = max(ROH_end),
    length    = ROH_end - ROH_start + 1,
    n_merged  = n(),  # how many original ROHs were merged
    .groups = "drop"
  )

# Save merged ROH table
write.csv(merged_roh, "~/Documents/SDZWA/alala_genomics/analysis/ROH/roh_combined_filtered_merged.csv")

# View a preview
head(merged_roh)



### recalculate FROH

# === 1. Provide chromosome lengths ===
chrom_lengths <- c(99062180,123451405,121534940,80112800,76278832,66143299,40978052,38659334,33574751,28816175,
                          23986915,23235100,22362767,23265108,21088201,17487702,16938386,15572589,12787533,13899114,
                          13093161,8696975,7971722,11344391,7312786,46055897,6387999,7803721,3496993,22021689,4714721,
                          736907,477022,104925,1321807,1784968,336161,251864,167149,337809,123330)

# === 2. Compute genome size per individual ===
genome_size <- sum(chrom_lengths)

# === 3. Compute total ROH length per individual for each threshold ===
thresholds <- c(1e6, 5e6, 10e6)

froh_list <- lapply(thresholds, function(thresh) {
  roh %>%
    filter(length >= thresh) %>%
    group_by(individual_ID) %>%
    summarise(total_roh_bp = sum(length), .groups = "drop") %>%
    mutate(
      threshold = paste0("FROH_", thresh / 1e6, "Mb"),
      FROH = total_roh_bp / genome_size
    )
}) %>%
  bind_rows()

froh_df <- froh_list %>%
  select(individual_ID, threshold, FROH) %>%
  pivot_wider(names_from = threshold, values_from = FROH)


# === 4. Save results ===
write.csv(froh_df, "~/Documents/SDZWA/alala_genomics/analysis/ROH/merged_FROH.csv")







