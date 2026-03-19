# Load libraries
library(ggplot2)
library(dplyr)
library(survival)
library(glmmTMB)
library(patchwork)  # for multi-panel layout
library(rstatix)


#---------------------------------
# Load data
#---------------------------------
ROH_data <- read.csv("~/Documents/SDZWA/alala_genomics/analysis/ROH/data/alala_medata.csv")

main_col <- "chartreuse4"

#---------------------------------
# Panel A: FROH vs hatch year
#---------------------------------
lm2 <- lm(FROH_1Mb ~ hatch_year, data = ROH_data)
summary(lm2)
p_hatch <- ggplot(ROH_data, aes(x = hatch_year, y = FROH_1Mb)) +
  geom_point(color = main_col, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  labs(x = "Hatch year", y = expression(F[ROH]), tag = "A") +
  coord_cartesian(ylim = c(0.1, 0.5)) +
  theme_bw(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18, hjust = -0.1))

#---------------------------------
# Panel B: Fped vs FROH
#---------------------------------
lm1 <- lm(FROH_1Mb ~ Fped, data = ROH_data)
summary(lm1)
p_fped <- ggplot(ROH_data, aes(x = Fped, y = FROH_1Mb)) +
  geom_point(color = main_col, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1.2) +
  coord_cartesian(xlim = c(0, 0.26), ylim = c(0.1, 0.5)) +
  labs(x = expression(F[PED]), y = expression(F[ROH]), tag = "B") +
  theme_bw(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18, hjust = -0.1))

#---------------------------------
# Panel C: Violin — embryo vs adult
#---------------------------------
alala_after2014 <- ROH_data %>%
  filter(hatch_year > 2013) %>%
  mutate(embryo = ifelse(lifespan > 0, "Hatched (n=31)", "Hatch failure (n=78)"))

test_embryo <- alala_after2014 %>%
  wilcox_test(FROH_1Mb ~ embryo) %>%
  add_significance()

test_embryo

p_violin_embryo <- ggplot(alala_after2014, aes(x = embryo, y = FROH_1Mb)) +
  geom_violin(fill = main_col, alpha = 1) +
  geom_jitter(width = 0.03, size = 1.5, alpha = 0.8, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 25, color = "black") +
  labs(x = "", y = expression(F[ROH]), tag = "C") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold", size = 18, hjust = -0.1)
  )

#---------------------------------
# Panel D: Survival residuals vs FROH
#---------------------------------
ROH_data_adults <- ROH_data %>% filter(lifespan > 0)


# 1. Fit Cox model
cox_model <- coxph(Surv(lifespan, status) ~ FROH_1Mb,
                           data = ROH_data_adults)
summary(cox_model)

# 2. Define representative FROH values (10th, 50th, 90th percentiles)
froh_vals <- quantile(ROH_data_adults$FROH_1Mb, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

newdata <- data.frame(FROH_1Mb = as.numeric(froh_vals))
rownames(newdata) <- c("0.23", "0.31", "0.38")

# 3. Generate survival curves
fit <- survfit(cox_model, newdata = newdata)

# 4. Convert survfit object to tidy dataframe
surv_df <- data.frame(
  time = rep(fit$time, ncol(fit$surv)),
  surv = as.vector(fit$surv),
  lower = as.vector(fit$lower),
  upper = as.vector(fit$upper),
  strata = rep(rownames(newdata), each = length(fit$time))
)

# 5. Clean up labels
surv_df$strata <- factor(surv_df$strata,
                         levels = c("0.23", "0.31", "0.38"))

# 6. Plot
p_survival_curves <- ggplot(surv_df, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  labs(
    x = "Time (days)",
    y = "Survival probability",
    color = expression(F[ROH]),
    fill = expression(F[ROH]),
    tag = "D"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = c(0.03, 0.03),
    legend.justification = c(0, 0),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.tag = element_text(face = "bold", size = 18, hjust = -0.1)
  )

p_survival_curves

#---------------------------------
# Panel E: Offspring per year vs FROH
#---------------------------------
m_rate <- glmmTMB(
  Total.no..of.nestlings ~ FROH_1Mb + (1|hatch_year),
  offset = log(lifespan / 365),
  family = nbinom2,
  data = ROH_data
)
summary(m_rate)
b <- fixef(m_rate)$cond

ROH_data <- ROH_data %>%
  mutate(rate = Total.no..of.nestlings / (lifespan / 365),
         pred_rate = exp(b[1] + b[2] * FROH_1Mb))

p_offspring <- ggplot(ROH_data, aes(x = FROH_1Mb, y = rate)) +
  geom_point(color = main_col, size = 2) +
  geom_line(aes(y = pred_rate), color = "black", linewidth = 1.2) +
  labs(x = expression(F[ROH]), y = "Offspring per year", tag = "E") +
  theme_bw(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18, hjust = -0.1))



#---------------------------------
# Panel F: Violin — breeder vs non-breeder
#---------------------------------
alala_age10 <- ROH_data %>%
  filter(lifespan > 10*365 & SampleType == "Bird") %>%
  mutate(breeder = ifelse(Total.no..of.nestlings > 0, "Breeder (n=51)", "Non-breeder (n=14)"))

test_breeder <- alala_age10 %>%
  wilcox_test(FROH_1Mb ~ breeder) %>%
  add_significance()

test_breeder

p_violin_breeder <- ggplot(alala_age10, aes(x = breeder, y = FROH_1Mb)) +
  geom_violin(fill = main_col, alpha = 1) +
  geom_jitter(width = 0.03, size = 1.5, alpha = 0.8, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 25, color = "black") +
  labs(x = "", y = expression(F[ROH]), tag = "F") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    plot.tag = element_text(face = "bold", size = 18, hjust = -0.1)
  )
#---------------------------------
# Combine into a 2x3 grid and save
#---------------------------------
final_fig <- (p_hatch | p_fped) / 
  (p_violin_embryo  | p_survival_curves) /
  (p_offspring | p_violin_breeder) +
  plot_annotation(tag_levels = NULL)

ggsave("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/FROH_combined_1Mb.pdf",
       final_fig, width = 8, height = 10)



mean(ROH_data$FROH_1Mb[ROH_data$SampleType=='Egg'])
mean(ROH_data$FROH_1Mb[ROH_data$SampleType=='Bird'])


final_fig <- (p_hatch | p_fped |p_violin_embryo) /
  (p_survival_curves | p_offspring | p_violin_breeder) +
  plot_annotation(tag_levels = NULL)

ggsave("~/Documents/SDZWA/alala_genomics/analysis/ROH/Plots/FROH_combined_1Mb_horiz.pdf",
       final_fig, width = 12, height = 7)

