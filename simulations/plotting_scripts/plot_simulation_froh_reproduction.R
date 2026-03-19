
library(glmmTMB)


dat_5e9 <- read.table("~/Documents/SDZWA/alala_genomics/simulations/data/sim_output_mu5e9_2.txt")
froh <- as.numeric(dat_5e9[1,])
offspring <- as.numeric(dat_5e9[2,]) 
age <- as.numeric(dat_5e9[3,])
df_5e9 <- data.frame(froh, offspring, age)

### negative binomial model with hatch year as random effect and lifespan as fixed effect
m_rate_5e9 <- glmmTMB(
  offspring ~ froh,
  offset = log(age),
  family = nbinom2, 
  data = df_5e9
)

summary(m_rate_5e9)
b <- fixef(m_rate_5e9)$cond

df_5e9 <- df_5e9 %>%
  mutate(rate = offspring / age,
         pred_rate = exp(b[1] + b[2] * froh))

plot_5e9 <- ggplot(df_5e9, aes(x = froh, y = rate)) +
  geom_point(color = 'chartreuse4', size = 2) +
  geom_line(aes(y = pred_rate), color = "black", linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 1.8)) +
  labs(x = expression(F[ROH]), y = "Offspring per year", tag = "") +
  theme_bw(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18, hjust = -0.1))

plot_5e9




dat_15e8 <- read.table("~/Documents/SDZWA/alala_genomics/simulations/data/sim_output_mu15e8_2.txt")
froh <- as.numeric(dat_15e8[1,])
offspring <- as.numeric(dat_15e8[2,]) 
age <- as.numeric(dat_15e8[3,])
df_15e8 <- data.frame(froh, offspring, age)

### negative binomial model with hatch year as random effect and lifespan as fixed effect
m_rate_15e8 <- glmmTMB(
  offspring ~ froh,
  offset = log(age),
  family = nbinom2, 
  data = df_15e8
)

summary(m_rate_15e8)
b <- fixef(m_rate_15e8)$cond

df_15e8 <- df_15e8 %>%
  mutate(rate = offspring / age,
         pred_rate = exp(b[1] + b[2] * froh))

plot_15e8 <- ggplot(df_15e8, aes(x = froh, y = rate)) +
  geom_point(color = 'chartreuse4', size = 2) +
  geom_line(aes(y = pred_rate), color = "black", linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 1.8)) +
  labs(x = expression(F[ROH]), y = "Offspring per year", tag = "") +
  theme_bw(base_size = 14) +
  theme(plot.tag = element_text(face = "bold", size = 18, hjust = -0.1))
plot_15e8





final_fig <- (plot_5e9 | plot_15e8) 

ggsave("~/Documents/SDZWA/alala_genomics/simulations/plots/sim_output.pdf",
       final_fig, width = 7, height = 3.5)
















