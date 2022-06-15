# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(insight)
library(bayestestR)
library(bayesplot)
library(cowplot)
library(brms)
library(broom)
library(ggpubr)
library(parameters)
library(multipanelfigure)
library(pvaluefunctions)
library(sanon)
#library(rstanarm)

# Load trial data
dat <- read.csv("/Users/blette/Downloads/outcomes.csv")

#######################################################################
# Replicate reported results
# Primary analysis
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup)
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup, correct = F)

dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
dat_primary$AgeFactor <- as.factor(dat_primary$AgeGroup)
mod_primary_rr <- glm(PrimaryEndpoint ~ Trt + AgeFactor,
                      data = dat_primary,
                      family = binomial(link = "log"))
mod_primary_rd <- glm(PrimaryEndpoint ~ Trt + AgeFactor,
                      data = dat_primary,
                      family = binomial(link = "identity"))

# Secondary analyses
dat_secondary1 <- dat[dat$SurviveM12 != 97, ]
mantelhaen.test(dat_secondary1$TreatRand, dat_secondary1$SurviveM12,
                dat_secondary1$AgeGroup, correct = F)

dat_secondary1$Trt <- abs(2 - dat_secondary1$TreatRand)
dat_secondary1$AgeFactor <- as.factor(dat_secondary1$AgeGroup)
mod_secondary_rr <- glm(SurviveM12 ~ Trt + AgeFactor,
                        data = dat_secondary1,
                        family = binomial(link = "log"))
mod_secondary_rd <- glm(SurviveM12 ~ Trt + AgeFactor,
                        data = dat_secondary1,
                        family = binomial(link = "identity"))

dat_secondary2 <- dat[dat$SurviveM12 != 97 & !is.na(dat$DeltaVabs), ]
summary(sanon(DeltaVabs ~ grp(TreatRand) + strt(AgeGroup),
              data = dat_secondary2))

# Make p-value function figures
# Primary outcome RR scale
pvf1 <- conf_dist(estimate = coef(mod_primary_rr)[2],
                  stderr = summary(mod_primary_rr)$coefficients[2, 2],
                  type = "logreg", plot_type = "p_val",
                  n_values = 1e4L,
                  null_values = c(0), trans = "exp", alternative = "two_sided",
                  log_yaxis = TRUE, cut_logyaxis = 0.05,
                  xlab = "Relative Benefit", ylab = "P-value (two-sided)",
                  ylab_sec = "P-value (one-sided)",
                  together = FALSE, plot_p_limit = 1 - 0.999,
                  plot_counternull = FALSE, plot = TRUE)

# Primary outcome RD scale
pvf2 <- conf_dist(estimate = 100*coef(mod_primary_rd)[2],
                  stderr = 100*summary(mod_primary_rd)$coefficients[2, 2],
                  type = "logreg", plot_type = "p_val",
                  n_values = 1e4L,
                  null_values = c(0), trans = "identity",
                  alternative = "two_sided",
                  log_yaxis = TRUE, cut_logyaxis = 0.05,
                  xlab = "Absolute Benefit Difference (%)",
                  ylab = "P-value (two-sided)",
                  ylab_sec = "P-value (one-sided)",
                  together = FALSE, plot_p_limit = 1 - 0.999,
                  plot_counternull = FALSE, plot = TRUE)

# Secondary outcome RR scale
pvf3 <- conf_dist(estimate = coef(mod_secondary_rr)[2],
                  stderr = summary(mod_secondary_rr)$coefficients[2, 2],
                  type = "logreg", plot_type = "p_val",
                  n_values = 1e4L,
                  null_values = c(0), trans = "exp", alternative = "two_sided",
                  log_yaxis = TRUE, cut_logyaxis = 0.05,
                  xlab = "Relative Benefit", ylab = "P-value (two-sided)",
                  ylab_sec = "P-value (one-sided)",
                  together = FALSE, plot_p_limit = 1 - 0.999,
                  plot_counternull = FALSE, plot = TRUE)

# Primary outcome RD scale
pvf4 <- conf_dist(estimate = 100*coef(mod_secondary_rd)[2],
                  stderr = 100*summary(mod_secondary_rd)$coefficients[2, 2],
                  type = "logreg", plot_type = "p_val",
                  n_values = 1e4L,
                  null_values = c(0), trans = "identity",
                  alternative = "two_sided",
                  log_yaxis = TRUE, cut_logyaxis = 0.05,
                  xlab = "Absolute Benefit Difference (%)",
                  ylab = "P-value (two-sided)",
                  ylab_sec = "P-value (one-sided)",
                  together = FALSE, plot_p_limit = 1 - 0.999,
                  plot_counternull = FALSE, plot = TRUE)


#######################################################################
# Bayesian re-analysis of primary outcome using flat priors

# SD = 100 ensures that the prior is uninformative
flat_prior <- prior(normal(0, 1000), class = "b")

thapca_flat <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
                   prior = flat_prior, family = "bernoulli", seed = 1234,
                   iter = 10000, chains = 8)

pred1 <- posterior_epred(thapca_flat,
                         newdata = mutate(thapca_flat$data, Trt = 1))

pred0 <- posterior_epred(thapca_flat,
                         newdata = mutate(thapca_flat$data, Trt = 0))

# First do relative risk scale
ratio <- apply(pred1, 1, mean) / apply(pred0, 1, mean)

# Some summaries used for plot legend
any_harm <- sum(ratio < 1) / length(ratio)
severe_harm <- sum(ratio < 1/1.25) / length(ratio)
rope <- sum(ratio < 1.05 & ratio > 1/1.05) / length(ratio)

# Some summaries used for separate table (see neuro-outcome-benefit-table.R)
fl_med_r <- median(ratio)
fl_low_r <- quantile(ratio, 0.025)
fl_up_r <- quantile(ratio, 0.975)
fl1.00 <- mean(ratio > 1)
fl1.10 <- mean(ratio > 1.1)
fl1.25 <- mean(ratio > 1.25)
fl1.50 <- mean(ratio > 1.5)
fl2.00 <- mean(ratio > 2)

# First figure panel
temp_plot <- ggplot(data.frame(ratio), aes(x = ratio)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(expand = c(0, 0), trans = "log",
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 6.2),
                  ylim = c(0, 1)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
  #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
  #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
  # 5), ],
  #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > 1/1.05 & x < 1.05)[
                   #c(1, 3, 6, 8), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 2.7, y = 0.8, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2))) +
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.905, ymax = 0.965,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.813, ymax = 0.873,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 2.2, xmax = 2.4, ymin = 0.813, ymax = 0.873,
           fill = "#3182BD") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.721, ymax = 0.781,
           fill = "#3182BD") +
  #annotate("segment", x = 2.42, xend = 2.42, y = 0.689, yend = 0.629) +
  #annotate("segment", x = 2.46, xend = 2.46, y = 0.689, yend = 0.629) +
  #annotate("segment", x = 2.50, xend = 2.50, y = 0.689, yend = 0.629) +
  #annotate("segment", x = 2.54, xend = 2.54, y = 0.689, yend = 0.629) +
  #annotate("segment", x = 2.58, xend = 2.58, y = 0.689, yend = 0.629) +
  geom_vline(xintercept = 1, color = "black", linetype = 1)

# Version that looks better when not using the p-value function panels
fig1_alt <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(expand = c(0, 0), trans = "log",
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 6.2),
                  ylim = c(0, 1)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
  #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
  #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
  # 5), ],
  #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > 1/1.05 & x < 1.05)[
                   #c(1, 3, 6, 8), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 2.7, y = 0.82, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2))) +
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.843, ymax = 0.873,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.803, ymax = 0.833,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 2.2, xmax = 2.4, ymin = 0.803, ymax = 0.833,
           fill = "#3182BD") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 0.763, ymax = 0.793,
           fill = "#3182BD") +
  #annotate("segment", x = 2.42, xend = 2.42, y = 0.723, yend = 0.753) +
  #annotate("segment", x = 2.46, xend = 2.46, y = 0.723, yend = 0.753) +
  #annotate("segment", x = 2.50, xend = 2.50, y = 0.723, yend = 0.753) +
  #annotate("segment", x = 2.54, xend = 2.54, y = 0.723, yend = 0.753) +
  #annotate("segment", x = 2.58, xend = 2.58, y = 0.723, yend = 0.753) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)

#######################################################################
# Primary outcome on RD scale
diff <- 100*(apply(pred1, 1, mean) - apply(pred0, 1, mean))

# Some summaries used for plot legend
any_harm <- sum(diff < 0) / length(diff)
severe_harm <- sum(diff < -5) / length(diff)
rope <- sum(diff < 1 & diff > -1) / length(diff)

# Some summaries used for separate table (see neuro-outcome-benefit-table.R)
fl_med_d <- median(diff)
fl_low_d <- quantile(diff, 0.025)
fl_up_d <- quantile(diff, 0.975)
fl0 <- mean(diff > 0)
fl2 <- mean(diff > 2)
fl5 <- mean(diff > 5)
fl10 <- mean(diff > 10)
fl15 <- mean(diff > 15)

# Second figure panel
temp_plot <- ggplot(data.frame(diff), aes(x = diff)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 0] <- 1
plot_data$barriers[plot_data$x < -5] <- 2

fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-10, 20, by = 5),
                     breaks = seq(-10, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-13, 33),
                  ylim = c(0, 0.095)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black", linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > -1 & x < 1)[
                   #c(1, 7, 21, 27), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 16.5, y = 0.075, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01")) +
                          #"\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0855, ymax = 0.0915,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0765, ymax = 0.0825,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 13.0, xmax = 14.5, ymin = 0.0765, ymax = 0.0825,
           fill = "#3182BD") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0675, ymax = 0.0735,
           fill = "#3182BD") +
  #annotate("segment", x = 14.65, xend = 14.65, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 14.95, xend = 14.95, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.25, xend = 15.25, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.55, xend = 15.55, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.85, xend = 15.85, y = 0.0585, yend = 0.0645) +
  geom_vline(xintercept = 0, color = "black", linetype = 1)

# Version that looks better when not using the p-value function panels
fig2_alt <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-10, 20, by = 5),
                     breaks = seq(-10, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-13, 33),
                  ylim = c(0, 0.095)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > -1 & x < 1)[
                   #c(1, 7, 21, 27), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 16.5, y = 0.077, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01")) +
                          #"\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0793, ymax = 0.0823,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0753, ymax = 0.0783,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 13.0, xmax = 14.5, ymin = 0.0753, ymax = 0.0783,
           fill = "#3182BD") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0713, ymax = 0.0743,
           fill = "#3182BD") +
  #annotate("segment", x = 14.65, xend = 14.65, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 14.95, xend = 14.95, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.25, xend = 15.25, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.55, xend = 15.55, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.85, xend = 15.85, y = 0.0673, yend = 0.0703) +
  geom_vline(xintercept = 0, color = "black", linetype = 1)

#######################################################################
# Bayesian re-analysis of secondary outcome

# Flat prior
thapca_flat2 <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
                    prior = flat_prior, family = "bernoulli",
                    seed = 1234, iter = 10000, chains = 8)

pred1 <- posterior_epred(thapca_flat2,
                         newdata = mutate(thapca_flat2$data, Trt = 1))

pred0 <- posterior_epred(thapca_flat2,
                         newdata = mutate(thapca_flat2$data, Trt = 0))

ratio <- apply(pred1, 1, mean) / apply(pred0, 1, mean)

# Some summaries used for plot legend
any_harm <- sum(ratio < 1) / length(ratio)
severe_harm <- sum(ratio < 1/1.25) / length(ratio)
rope <- sum(ratio < 1.05 & ratio > 1/1.05) / length(ratio)

# Some summaries used for separate table (see survival-outcome-benefit-table.R)
fl_med_r <- median(ratio)
fl_low_r <- quantile(ratio, 0.025)
fl_up_r <- quantile(ratio, 0.975)
fl1.00 <- mean(ratio > 1)
fl1.10 <- mean(ratio > 1.1)
fl1.25 <- mean(ratio > 1.25)
fl1.50 <- mean(ratio > 1.5)
fl2.00 <- mean(ratio > 2)

# Third figure panel
temp_plot <- ggplot(data.frame(ratio), aes(x = ratio)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig3 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(expand = c(0, 0), trans = "log",
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 6.2),
                  ylim = c(0, 2)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > 1/1.05 & x < 1.05)[
                   #c(1, 7, 19, 25), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 2.7, y = 1.6, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01")) +
                          #"\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.81, ymax = 1.93,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.626, ymax = 1.746,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 2.2, xmax = 2.4, ymin = 1.626, ymax = 1.746,
           fill = "#3182BD") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.442, ymax = 1.562,
           fill = "#3182BD") +
  #annotate("segment", x = 2.42, xend = 2.42, y = 1.378, yend = 1.258) +
  #annotate("segment", x = 2.46, xend = 2.46, y = 1.378, yend = 1.258) +
  #annotate("segment", x = 2.50, xend = 2.50, y = 1.378, yend = 1.258) +
  #annotate("segment", x = 2.54, xend = 2.54, y = 1.378, yend = 1.258) +
  #annotate("segment", x = 2.58, xend = 2.58, y = 1.378, yend = 1.258) +
  geom_vline(xintercept = 1, color = "black", linetype = 1)

# Version that looks better when not using the p-value function panels
fig3_alt <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(expand = c(0, 0), trans = "log",
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 6.2),
                  ylim = c(0, 2)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > 1/1.05 & x < 1.05)[
                   #c(1, 7, 19, 25), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 2.7, y = 1.637, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01")) +
                          #"\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.686, ymax = 1.746,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.606, ymax = 1.666,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 2.2, xmax = 2.4, ymin = 1.606, ymax = 1.666,
           fill = "#3182BD") +
  annotate("rect", xmin = 2.4, xmax = 2.6, ymin = 1.526, ymax = 1.586,
           fill = "#3182BD") +
  #annotate("segment", x = 2.42, xend = 2.42, y = 1.446, yend = 1.506) +
  #annotate("segment", x = 2.46, xend = 2.46, y = 1.446, yend = 1.506) +
  #annotate("segment", x = 2.50, xend = 2.50, y = 1.446, yend = 1.506) +
  #annotate("segment", x = 2.54, xend = 2.54, y = 1.446, yend = 1.506) +
  #annotate("segment", x = 2.58, xend = 2.58, y = 1.446, yend = 1.506) +
  geom_vline(xintercept = 1, color = "black", linetype = 1)

#######################################################################
# Secondary outcome on RD scale
diff <- 100*(apply(pred1, 1, mean) - apply(pred0, 1, mean))

# Some summaries used for plot legend
any_harm <- sum(diff < 0) / length(diff)
severe_harm <- sum(diff < -5) / length(diff)
rope <- sum(diff < 1 & diff > -1) / length(diff)

# Some summaries used for separate table (see survival-outcome-benefit-table.R)
fl_med_d <- median(diff)
fl_low_d <- quantile(diff, 0.025)
fl_up_d <- quantile(diff, 0.975)
fl0 <- mean(diff > 0)
fl2 <- mean(diff > 2)
fl5 <- mean(diff > 5)
fl10 <- mean(diff > 10)
fl15 <- mean(diff > 15)

# Fourth figure panel
temp_plot <- ggplot(data.frame(diff), aes(x = diff)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 0] <- 1
plot_data$barriers[plot_data$x < -5] <- 2

fig4 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-10, 20, by = 5),
                     breaks = seq(-10, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-13, 33),
                  ylim = c(0, 0.095)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > -1 & x < 1)[
                   #c(1, 6, 18, 23), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 16.5, y = 0.075, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2))) +
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0855, ymax = 0.0915,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0765, ymax = 0.0825,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 13.0, xmax = 14.5, ymin = 0.0765, ymax = 0.0825,
           fill = "#3182BD") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0675, ymax = 0.0735,
           fill = "#3182BD") +
  #annotate("segment", x = 14.65, xend = 14.65, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 14.95, xend = 14.95, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.25, xend = 15.25, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.55, xend = 15.55, y = 0.0585, yend = 0.0645) +
  #annotate("segment", x = 15.85, xend = 15.85, y = 0.0585, yend = 0.0645) +
  geom_vline(xintercept = 0, color = "black", linetype = 1)

# Version that looks better when not using the p-value function panels
fig4_alt <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-10, 20, by = 5),
                     breaks = seq(-10, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-13, 33),
                  ylim = c(0, 0.095)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  #geom_segment(inherit.aes = FALSE,
               #data =
                 #subset(plot_data, x > -1 & x < 1)[
                   #c(1, 6, 18, 23), ],
               #aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 16.5, y = 0.077, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2))) +
                          #round(severe_harm, 2),
                          #"\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0793, ymax = 0.0823,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0753, ymax = 0.0783,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 13.0, xmax = 14.5, ymin = 0.0753, ymax = 0.0783,
           fill = "#3182BD") +
  annotate("rect", xmin = 14.5, xmax = 16.0, ymin = 0.0713, ymax = 0.0743,
           fill = "#3182BD") +
  #annotate("segment", x = 14.65, xend = 14.65, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 14.95, xend = 14.95, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.25, xend = 15.25, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.55, xend = 15.55, y = 0.0673, yend = 0.0703) +
  #annotate("segment", x = 15.85, xend = 15.85, y = 0.0673, yend = 0.0703) +
  geom_vline(xintercept = 0, color = "black", linetype = 1)


############################################################################
# Combine figures

# First do Bayesian analyses on RR scale (first version of fig)
# Both figures exported as landscape PDFs with 6 in x 10 in dimensions
figure <- multi_panel_figure(columns = 2, rows = 1, width = 240,
                             height = 125)
figure <- fill_panel(figure, fig1_alt)
figure <- fill_panel(figure, fig3_alt)

# Output pdf of RR figure, dims may need to be changed
pdf("flat-prior-RR-figure.pdf", width = 10, height = 6)
figure
dev.off()

# Then Bayesian analyses on RD scale (first version of fig)
figure2 <- multi_panel_figure(columns = 2, rows = 1, width = 250,
                              height = 125)
figure2 <- fill_panel(figure2, fig2_alt)
figure2 <- fill_panel(figure2, fig4_alt)

# Output pdf of RD figure, dims may need to be changed
pdf("flat-prior-RD-figure.pdf", width = 10, height = 6)
figure2
dev.off()

# All RR scale analyses (second version of fig)
figure3 <- multi_panel_figure(columns = 2, rows = 2, width = 240,
                              height = 125)
figure3 <- fill_panel(figure3,
                      pvf1$plot +
                        geom_hline(yintercept = c(0.05, 0.1, 0.2),
                                   linetype = 2) +
                        theme(axis.text.x = element_text(size = 8),
                              axis.text.y = element_text(size = 8),
                              axis.title.y.left = element_text(size = 10),
                              axis.title.y.right = element_text(size = 10),
                              axis.title.x = element_text(size = 10)))
figure3 <- fill_panel(figure3, fig1)
figure3 <- fill_panel(figure3,
                      pvf3$plot +
                        geom_hline(yintercept = c(0.05, 0.1, 0.2),
                                   linetype = 2) +
                        theme(axis.text.x = element_text(size = 8),
                              axis.text.y = element_text(size = 8),
                              axis.title.y.left = element_text(size = 10),
                              axis.title.y.right = element_text(size = 10),
                              axis.title.x = element_text(size = 10)))
figure3 <- fill_panel(figure3, fig3)

# Output pdf of RR figure, dims may need to be changed
pdf("pvf-flat-prior-RR-figure.pdf", width = 9.75, height = 6)
figure3
dev.off()

# All RD scale analyses (second version of fig)
figure4 <- multi_panel_figure(columns = 2, rows = 2, width = 240,
                              height = 125)
figure4 <- fill_panel(figure4,
                      pvf2$plot +
                        geom_hline(yintercept = c(0.05, 0.1, 0.2),
                                   linetype = 2) +
                        theme(axis.text.x = element_text(size = 8),
                              axis.text.y = element_text(size = 8),
                              axis.title.y.left = element_text(size = 10),
                              axis.title.y.right = element_text(size = 10),
                              axis.title.x = element_text(size = 10)))
figure4 <- fill_panel(figure4, fig2)
figure4 <- fill_panel(figure4,
                      pvf4$plot +
                        geom_hline(yintercept = c(0.05, 0.1, 0.2),
                                   linetype = 2) +
                        theme(axis.text.x = element_text(size = 8),
                              axis.text.y = element_text(size = 8),
                              axis.title.y.left = element_text(size = 10),
                              axis.title.y.right = element_text(size = 10),
                              axis.title.x = element_text(size = 10)))
figure4 <- fill_panel(figure4, fig4)

# Output pdf of RR figure, dims may need to be changed
pdf("pvf-flat-prior-RD-figure.pdf", width = 9.75, height = 6)
figure4
dev.off()


##########################################################################
# Run model diagnostics

# Borrow diagnostics function from COVID Steroid 2 Trial analysis
# MCMC diagnostics
diag_fit <- function(fit, extra_groups = NULL, bars = TRUE) {
  
  # Print fit (includes Rhats and bulk and tail ESS)
  print(fit)
  if (any(rhat(fit) > 1.01)) {
    warning("One or more Rhats > 1.01")
  } else {
    message("All Rhats <= 1.01")
  }
  
  nam <- names(fit$fit)
  walk(nam, ~{print(mcmc_trace(fit, pars = .x))})
  walk(nam, ~{print(mcmc_dens_overlay(fit, pars = .x))})
  if (bars) {
    print(pp_check(fit, nsamples = 100, type = "bars"))
    print(pp_check(fit, nsamples = 100, type = "bars_grouped", group = "Trt"))
  } else {
    print(pp_check(fit, nsamples = 100, type = "dens_overlay"))
    print(pp_check(fit, nsamples = 100, type = "dens_overlay_grouped",
                   group = "Trt"))
  }
  print(pp_check(fit, type = "stat_grouped", stat = "mean", group = "Trt"))
  
  if (!is.null(extra_groups)) {
    
    for (i in seq_along(extra_groups)) {
      
      if (bars) {
        print(pp_check(fit, nsamples = 100, type = "bars_grouped",
                       group = extra_groups[i]))
      } else {
        print(pp_check(fit, nsamples = 100, type = "dens_overlay_grouped",
                       group = extra_groups[i]))
      }
      print(pp_check(fit, type = "stat_grouped", stat = "mean",
                     group = extra_groups[i]))
      
    }
    
  }
  
  print(loo(fit, cores = 4, reloo = TRUE)) # RELOO IF NECESSARY
  invisible(fit)
  
}

diag_fit(thapca_flat)
diag_fit(thapca_flat2)
