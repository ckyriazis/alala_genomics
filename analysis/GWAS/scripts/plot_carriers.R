
library(dplyr)
library(ggplot2)
library(patchwork)
library(logistf)
library(tidyverse)

ROH_data <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/alala_medata.csv")
repro_data <-  read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/alala_repro_history_summary_2025.09.23.csv")
data_merged <- merge(ROH_data, repro_data, by="Sample")

ROH_data_adults <- ROH_data[ROH_data$lifespan>0 & ROH_data$hatch_year > 2012,]


library(patchwork)

dat <- ROH_data_adults %>%
  mutate(
    offspring_per_year = Total.no..of.nestlings / (lifespan / 365.25),
    DLG1_carrier = case_when(
      DLG1 == 0 ~ "Non-carrier",
      DLG1 %in% c(1, 2) ~ "Carrier",
      TRUE ~ NA_character_
    ),
    NEO1_carrier = case_when(
      NEO1 == 0 ~ "Non-carrier",
      NEO1 %in% c(1, 2) ~ "Carrier",
      TRUE ~ NA_character_
    )
  )


wilcox.test(offspring_per_year ~ DLG1_carrier, data = dat)
wilcox.test(offspring_per_year ~ NEO1_carrier, data = dat)


dat %>%
  group_by(DLG1_carrier) %>%
  summarise(mean_offspring_per_year = mean(offspring_per_year, na.rm = TRUE))

dat %>%
  group_by(NEO1_carrier) %>%
  summarise(mean_offspring_per_year = mean(offspring_per_year, na.rm = TRUE))


dat_subset <- dat %>% filter(!is.na(DLG1_carrier))
dim(dat_subset)

dat_subset <- dat %>% filter(!is.na(NEO1_carrier))
dim(dat_subset)

table(dat$DLG1)
table(dat$NEO1)

# Create labeled x variables with sample sizes
dlg1_labels <- dat %>%
  filter(!is.na(DLG1_carrier)) %>%
  count(DLG1_carrier) %>%
  mutate(label = paste0(DLG1_carrier, "\n(n=", n, ")"))

neo1_labels <- dat %>%
  filter(!is.na(NEO1_carrier)) %>%
  count(NEO1_carrier) %>%
  mutate(label = paste0(NEO1_carrier, "\n(n=", n, ")"))

dat <- dat %>%
  left_join(dlg1_labels %>% select(DLG1_carrier, DLG1_label = label), by = "DLG1_carrier") %>%
  left_join(neo1_labels %>% select(NEO1_carrier, NEO1_label = label), by = "NEO1_carrier")

p1 <- ggplot(dat %>% filter(!is.na(DLG1_carrier)), 
             aes(x = DLG1_label, y = offspring_per_year)) +
  geom_violin(fill = "orange3") +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6, fill = "gray40", color = "gray40") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  labs(x = "DLG1", y = "Offspring per year") +
  theme_classic()

p2 <- ggplot(dat %>% filter(!is.na(NEO1_carrier)), 
             aes(x = NEO1_label, y = offspring_per_year)) +
  geom_violin(fill = "orange3") +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.6, fill = "gray40", color = "gray40") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  labs(x = "NEO1", y = "Offspring per year") +
  theme_classic()

p1 + p2

ggsave("~/Documents/SDZWA/alala_genomics/analysis/GWAS/Plots/carriers_plot_2012.pdf", p1 + p2, width = 7, height = 3)




