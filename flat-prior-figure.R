# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(insight)
library(bayestestR)
library(cowplot)
library(brms)
library(broom)
library(ggpubr)
library(parameters)
library(sanon)
library(rstanarm)

# Load trial data
dat <- read.csv("../outcomes.csv")

#######################################################################
# Replicate reported results
# Primary analysis
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup)
mantelhaen.test(dat_primary$TreatRand, dat_primary$PrimaryEndpoint,
                dat_primary$AgeGroup, correct = F)

mod_primary <- glm(PrimaryEndpoint ~ TreatRand + factor(AgeGroup),
                   data = dat_primary,
                   family = binomial(link = "log"))

# Secondary analyses
dat_secondary1 <- dat[dat$SurviveM12 != 97, ]
mantelhaen.test(dat_secondary1$TreatRand, dat_secondary1$SurviveM12,
                dat_secondary1$AgeGroup, correct = F)

dat_secondary2 <- dat[dat$SurviveM12 != 97 & !is.na(dat$DeltaVabs), ]
summary(sanon(DeltaVabs ~ grp(TreatRand) + strt(AgeGroup),
              data = dat_secondary2))

#######################################################################
# Bayesian re-analysis of primary outcome

# Flat prior analysis
# SD = 100 ensures that the prior is uninformative
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
thapca_flat <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                        data = dat_primary,
                        prior_intercept = normal(0, 100),
                        prior = normal(location = c(0, 0),
                                       scale = c(100, 100)),
                        family = binomial(link = "log"), seed = 1234)

posteriors <- insight::get_parameters(thapca_flat, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(exp(posteriors$Trt) < 1) / nrow(posteriors)
severe_harm <- sum(exp(posteriors$Trt) < 1/1.25) / nrow(posteriors)
rope <- sum(exp(posteriors$Trt) < 1.05 & exp(posteriors$Trt) > 1/1.05) /
  nrow(posteriors)

# First figure panel
posteriors$exptrt <- exp(posteriors$Trt)
temp_plot <- ggplot(posteriors, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Risk", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(0.5, 4, by = 0.5),
                     breaks = seq(0.5, 4, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 4.5),
                  ylim = c(0, 1)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
  #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
  #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
  # 5), ],
  #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  geom_segment(inherit.aes = FALSE,
               data =
                 subset(plot_data, x > 1/1.05 & x < 1.05)[
                   c(1, 3, 6, 8), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 3.5, y = 0.9, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) = ",
                          round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.943, ymax = 0.973,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.903, ymax = 0.933,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 3.25, xmax = 3.35, ymin = 0.903, ymax = 0.933,
           fill = "#3182BD") +
  annotate("rect", xmin = 3.35, xmax = 3.45, ymin = 0.863, ymax = 0.893,
           fill = "#3182BD") +
  annotate("segment", x = 3.36, xend = 3.36, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.38, xend = 3.38, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.40, xend = 3.40, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.42, xend = 3.42, y = 0.853, yend = 0.823) +
  annotate("segment", x = 3.44, xend = 3.44, y = 0.853, yend = 0.823) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)

#######################################################################
# Bayesian re-analysis of secondary outcome

# Flat prior
dat_secondary1$Trt <- abs(2 - dat_secondary1$TreatRand)
thapca_flat2 <- stan_glm(SurviveM12 ~ Trt + AgeGroup,
                         data = dat_secondary1,
                         prior_intercept = normal(0, 100),
                         prior = normal(location = c(0, 0),
                                        scale = c(100, 100)),
                         family = binomial(link = "log"), seed = 1234)

posteriors2 <- insight::get_parameters(thapca_flat2, iterations = 10^5)

# Some summaries used for plot legend
any_harm <- sum(exp(posteriors2$Trt) < 1) / nrow(posteriors2)
severe_harm <- sum(exp(posteriors2$Trt) < 1/1.25) / nrow(posteriors2)
rope <- sum(exp(posteriors2$Trt) < 1.05 &
            exp(posteriors2$Trt) > 1/1.05) / nrow(posteriors2)

# Second figure panel
posteriors2$exptrt <- exp(posteriors2$Trt)
temp_plot <- ggplot(posteriors2, aes(x = exptrt)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]
plot_data$barriers <- 0
plot_data$barriers[plot_data$x < 1] <- 1
plot_data$barriers[plot_data$x < 1/1.25] <- 2

fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers)) +
  geom_line() +
  labs(x = "Relative Risk", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(0.5, 2.5, by = 0.5),
                     breaks = seq(0.5, 2.5, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 2.5),
                  ylim = c(0, 2)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  #geom_segment(inherit.aes = FALSE,
  #data = subset(plot_data, x > 1/1.05 & x < 1.05)[
  #seq(1, nrow(subset(plot_data, x > 1/1.05 & x < 1.05)),
  # 5), ],
  #aes(x = x, y = 0, xend = x, yend = y)) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  scale_fill_brewer(guide = "none") +
  geom_segment(inherit.aes = FALSE,
               data =
                 subset(plot_data, x > 1/1.05 & x < 1.05)[
                   c(1, 7, 21, 28), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  annotate("text", x = 2.1, y = 1.6, size = 3, hjust = 0,
           label = paste0("P(Benefit) = ", 1 - round(any_harm, 2),
                          "\n", "P(Harm) = ", round(any_harm, 2),
                          "\n", "P(Severe Harm) < 0.01",
                          #round(severe_harm, 2),
                          "\n", "ROPE = ", round(rope, 2))) +
  annotate("rect", xmin = 2.03, xmax = 2.08, ymin = 1.68, ymax = 1.75,
           fill = "#DEEBF7") +
  annotate("rect", xmin = 2.03, xmax = 2.08, ymin = 1.6, ymax = 1.67,
           fill = "#9ECAE1") +
  annotate("rect", xmin = 1.98, xmax = 2.03, ymin = 1.6, ymax = 1.67,
           fill = "#3182BD") +
  annotate("rect", xmin = 2.03, xmax = 2.08, ymin = 1.52, ymax = 1.59,
           fill = "#3182BD") +
  annotate("segment", x = 2.035, xend = 2.035, y = 1.44, yend = 1.51) +
  annotate("segment", x = 2.045, xend = 2.045, y = 1.44, yend = 1.51) +
  annotate("segment", x = 2.055, xend = 2.055, y = 1.44, yend = 1.51) +
  annotate("segment", x = 2.065, xend = 2.065, y = 1.44, yend = 1.51) +
  annotate("segment", x = 2.075, xend = 2.075, y = 1.44, yend = 1.51) +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1)

ggarrange(fig1, fig2, ncol = 2)
