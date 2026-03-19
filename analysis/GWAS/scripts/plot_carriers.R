
library(dplyr)
library(ggplot2)
library(patchwork)

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
    HatchRateCond = Hatch_years /Egg_years,
    FROH_1Mb = mean(FROH_1Mb),
    FROH_10Mb = mean(FROH_10Mb),
    hatch_year = mean(hatch_year), 
    sex = Sex[1],
    
    carrier_DLG1 = sum(DLG1, na.rm = TRUE) / N_years,
    carrier_NEO1 = sum(NEO1, na.rm = TRUE) / N_years,
    carrier_all = sum(DLG1 + NEO1, na.rm = TRUE) / N_years
  )
indiv_summary
dim(indiv_summary)


indiv_summary_subset <- indiv_summary[indiv_summary$Opp_years>0 ,]
indiv_summary_subset_2012 <- indiv_summary[indiv_summary$Opp_years>0 & indiv_summary$hatch_year>2012 ,]

dim(indiv_summary_subset)


carrier2 <- indiv_summary_subset[indiv_summary_subset$carrier_all==2,]
carrier1 <- indiv_summary_subset[indiv_summary_subset$carrier_all==1,]
carrier0 <- indiv_summary_subset[indiv_summary_subset$carrier_all==0,]

### annual egg hatch rates by carrier status
mean(carrier2$HatchRate)
mean(carrier1$HatchRate)
mean(carrier0$HatchRate)

### conditional rates
mean(carrier2$HatchRate)/mean(carrier2$EggRate)
mean(carrier1$HatchRate)/mean(carrier1$EggRate)
mean(carrier0$HatchRate)/mean(carrier0$EggRate)



### plot only annual egg hatch rate 
plot <- ggplot(indiv_summary_subset_2012, 
             aes(x = factor(carrier_all), y = HatchRate, fill = factor(carrier_all))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(x = "Recessive lethal carrier status", y = "Average annual hatch rate") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 0.7))
plot


ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/Plots/carrier_hatch_rate_2012.pdf",
       plot, width = 4, height = 4)






### plot egg rate and conditional hatch rate for subset and all data

p1 <- ggplot(indiv_summary_subset, 
             aes(x = factor(carrier_all), y = EggRate, fill = factor(carrier_all))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(x = "Recessive lethal carrier status", y = "Egg laying rate") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("All individuals (n=84)")

p2 <- ggplot(indiv_summary_subset, 
               aes(x = factor(carrier_all), y = HatchRateCond, fill = factor(carrier_all))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(x = "Recessive lethal carrier status", y = "Conditional hatch rate") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("All individuals (n=84)")

p3 <- ggplot(indiv_summary_subset_2012, 
             aes(x = factor(carrier_all), y = EggRate, fill = factor(carrier_all))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(x = "Recessive lethal carrier status", y = "Egg laying rate") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1))+
  ggtitle("Post-2012 individuals (n=25)")

### plot only annual egg hatch rate 
p4 <- ggplot(indiv_summary_subset_2012, 
             aes(x = factor(carrier_all), y = HatchRateCond, fill = factor(carrier_all))) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
  labs(x = "Recessive lethal carrier status", y = "Conditional hatch rate") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(0, 1)) + 
  ggtitle("Post-2012 individuals (n=25)")


final_fig <- (p1 + p2) / (p3 + p4)
final_fig



ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/Plots/carrier_hatch_egg_rate_all.pdf",
       final_fig, width = 8, height = 7)


