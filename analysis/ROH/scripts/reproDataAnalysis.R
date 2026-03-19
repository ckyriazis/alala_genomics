
library(dplyr)

ROH_data <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/alala_medata.csv")
repro_data <-  read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/alala_repro_history_summary_2025.09.23.csv")
data_merged <- merge(ROH_data, repro_data, by="Sample")


indiv_summary <- data_merged %>%
  group_by(Sample) %>%
  summarise(
    # opportunity
    N_years = n(),                               # how many breeding seasons observed
    Opp_years = sum(MatesPerYr > 0, na.rm = TRUE),             # years given a mate
    
    # raw outputs
    Egg_years = sum(EggLaid, na.rm = TRUE),
    Fertile_years = sum(FertileEggConfirmed, na.rm = TRUE),
    Hatch_years = sum(EggHatched, na.rm = TRUE),
    
    # rates conditional on opportunity
    EggRate = Egg_years / Opp_years,
    FertileRate = Fertile_years / Opp_years,
    HatchRate = Hatch_years / Opp_years, 
    FROH_1Mb = mean(FROH_1Mb),
    FROH_10Mb = mean(FROH_10Mb),
    hatch_year = mean(hatch_year), 
    sex = Sex[1],
    
    carrier_NEO1 = sum(NEO1, na.rm = TRUE) / N_years,
    carrier_DLG1 = sum(DLG1, na.rm = TRUE) / N_years,
    carrier_all = sum(NEO1 + DLG1, na.rm = TRUE) / N_years
  )
indiv_summary


dim(indiv_summary[indiv_summary$Opp_years>0,])


pdf("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/alala_repro_results.pdf", width=6, height=5)
par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(indiv_summary$FROH_1Mb, indiv_summary$EggRate, col = "chartreuse4", pch=19, xlab = expression(' F'[ROH]), ylab = "Egg laying rate")
lm <- lm(EggRate~FROH_1Mb, data = indiv_summary)
summary(lm)
abline(a = lm$coefficients[1], b = lm$coefficients[2])


plot(indiv_summary$FROH_1Mb, indiv_summary$HatchRate, col = "chartreuse4", pch=19, xlab = expression(' F'[ROH]), ylab = "Egg hatch rate")
lm <- lm(HatchRate~FROH_1Mb, data = indiv_summary)
summary(lm)
abline(a = lm$coefficients[1], b = lm$coefficients[2])


plot(indiv_summary$FROH_10Mb, indiv_summary$EggRate, col = "chartreuse4", pch=19, xlab = expression(' F'[ROH]), ylab = "Egg laying rate")
lm <- lm(EggRate~FROH_10Mb, data = indiv_summary)
summary(lm)
abline(a = lm$coefficients[1], b = lm$coefficients[2])


plot(indiv_summary$FROH_10Mb, indiv_summary$HatchRate, col = "chartreuse4", pch=19, xlab = expression(' F'[ROH]), ylab = "Egg hatch rate")
lm <- lm(HatchRate~FROH_10Mb, data = indiv_summary)
summary(lm)
abline(a = lm$coefficients[1], b = lm$coefficients[2])


dev.off()




### Get minimum and maximum reproductive age 

# Step 1: Calculate age and filter for successful egg laying
breeding_data <- repro_data %>%
  mutate(Age = BreedingSeason - HatchYr) %>%
  filter(EggLaid == 1)

# Step 2: Get min/max breeding ages for each individual using base R aggregate
breeding_ages <- aggregate(
  Age ~ Sample + Studbook + Sex, 
  data = breeding_data, 
  FUN = function(x) c(MinBreedingAge = min(x), MaxBreedingAge = max(x))
)

# Unpack the matrix column into separate columns
breeding_ages <- data.frame(
  Sample = breeding_ages$Sample,
  Studbook = breeding_ages$Studbook,
  Sex = breeding_ages$Sex,
  MinBreedingAge = breeding_ages$Age[, "MinBreedingAge"],
  MaxBreedingAge = breeding_ages$Age[, "MaxBreedingAge"]
)

# Print individual results
print("Individual-level breeding ages:")
print(breeding_ages)


# Calculate summary statistics stratified by sex using base R
age_summary_by_sex <- do.call(rbind, lapply(split(breeding_ages, breeding_ages$Sex), function(x) {
  data.frame(
    Sex = unique(x$Sex),
    N_Individuals = nrow(x),
    Min_Age_Overall = min(x$MinBreedingAge),
    Max_Age_Overall = max(x$MaxBreedingAge),
    Mean_Min_Age = mean(x$MinBreedingAge),
    SD_Min_Age = sd(x$MinBreedingAge),
    Mean_Max_Age = mean(x$MaxBreedingAge),
    SD_Max_Age = sd(x$MaxBreedingAge)
  )
}))

print("\nSummary statistics by sex:")
print(age_summary_by_sex)

# Detailed summary
cat("\n=== Breeding Age Range by Sex ===\n")
for(sex in unique(age_summary_by_sex$Sex)) {
  sex_data <- age_summary_by_sex[age_summary_by_sex$Sex == sex, ]
  cat(sprintf("\n%s individuals (n=%d):\n", sex, sex_data$N_Individuals))
  cat(sprintf("  Minimum breeding age observed: %d years\n", sex_data$Min_Age_Overall))
  cat(sprintf("  Maximum breeding age observed: %d years\n", sex_data$Max_Age_Overall))
  cat(sprintf("  Average age at first successful breeding: %.1f ± %.1f years\n", 
              sex_data$Mean_Min_Age, sex_data$SD_Min_Age))
  cat(sprintf("  Average age at last successful breeding: %.1f ± %.1f years\n", 
              sex_data$Mean_Max_Age, sex_data$SD_Max_Age))
}



