library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(binom)


#------------------------------------------------------------
# LOAD DATA
#------------------------------------------------------------
df <- fread("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/alala_medata.csv")
df <- df %>% filter(grepl("AL", Sample))

# Parse dates
df <- df %>%
  mutate(
    Hatch_date = mdy(`Hatch date`),
    Death_date = mdy(`Death date; blank = alive`),
    # If Death_date is NA, assume still alive (use today's date or a fixed date)
    Death_date = if_else(is.na(Death_date), as.Date("2024-12-31"), Death_date)
  )

#------------------------------------------------------------
# CREATE TIME SERIES OF LIVING POPULATION
#------------------------------------------------------------

# Get range of years to analyze
year_range <- seq(min(df$hatch_year, na.rm = TRUE), 2024, by = 1)


# For each year, calculate which individuals were alive
allele_freq_living <- lapply(year_range, function(yr) {
  date_point <- as.Date(paste0(yr, "-12-31"))
  
  alive <- df %>%
    filter(
      !is.na(NEO1),
      !is.na(DLG1),
      Hatch_date <= date_point,
      Death_date >= date_point
    )
  
  n_ind <- nrow(alive)
  n_alleles <- 2 * n_ind  # Total number of alleles
  
  # Count derived alleles
  NEO1_count <- sum(alive$NEO1)
  DLG1_count <- sum(alive$DLG1)

  
  # Calculate exact binomial CIs
  NEO1_ci <- binom.confint(NEO1_count, n_alleles, method = "exact")
  DLG1_ci <- binom.confint(DLG1_count, n_alleles, method = "exact")

  data.frame(
    year = yr,
    n_ind = n_ind,
    NEO1_freq = NEO1_count / n_alleles,
    NEO1_lower = NEO1_ci$lower,
    NEO1_upper = NEO1_ci$upper,
    DLG1_freq = DLG1_count / n_alleles,
    DLG1_lower = DLG1_ci$lower,
    DLG1_upper = DLG1_ci$upper
  )
}) %>% bind_rows()

# Convert to long format
allele_freq_long <- allele_freq_living %>%
  pivot_longer(
    cols = c(NEO1_freq, DLG1_freq),
    names_to = "allele",
    values_to = "derived_freq"
  ) %>%
  mutate(
    lower = if_else(allele == "NEO1_freq", NEO1_lower, DLG1_lower),
    upper = if_else(allele == "NEO1_freq", NEO1_upper, DLG1_upper)
  )

# Plot with error bars
p <- ggplot(allele_freq_long, aes(
  x = year,
  y = derived_freq,
  color = allele,
  fill = allele
)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n_ind), alpha = 0.8) +
  scale_color_manual(
    values = c(NEO1_freq = "#009E73", DLG1_freq = "#D55E00"),
    labels = c("DLG1_freq" = "DLG1", "NEO1_freq" = "NEO1")
  ) +
  scale_fill_manual(
    values = c(NEO1_freq = "#009E73", DLG1_freq = "#D55E00"),
    guide = "none"
  ) +
  scale_size_continuous(name = "Sample size") +
  coord_cartesian(ylim = c(0, 0.4)) +
  labs(
    title = "",
    x = "Year",
    y = "Lethal haplotype frequency",
    color = "Lethal allele"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

print(p)


#------------------------------------------------------------
# SAVE
#------------------------------------------------------------
ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/plots/allele_frequency_living_population.png", 
       p, width = 8, height = 4, dpi = 300)

#------------------------------------------------------------
# OPTIONAL: Print summary statistics
#------------------------------------------------------------
cat("\nAllele frequency in current living population:\n")
current_living <- allele_freq_living %>% 
  filter(year == max(year))
print(current_living)

cat("\nChange in allele frequency over time:\n")
change_summary <- allele_freq_living %>%
  filter(year == min(year) | year == max(year)) %>%
  select(year, NEO1_freq, DLG1_freq)
print(change_summary)











