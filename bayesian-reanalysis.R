# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(insight)
library(bayestestR)
library(cowplot)
library(brms)
library(broom)
library(parameters)
library(sanon)

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
# SD = 1000 ensures that the prior is uninformative
flat_prior <- prior(normal(0, 1000), class = "b")

# Fit a Bayesian GLM (how best to handle stratified randomization?)
#art_flat <- brm(PrimaryEndpoint ~ TreatRand, data = dat_primary,
                #prior = flat_prior,
                #family = "bernoulli", seed = 123)
thapca_flat <- brm(PrimaryEndpoint ~ TreatRand + AgeGroup,
                data = dat_primary,
                prior = flat_prior,
                family = "bernoulli", seed = 123)
posteriors <- insight::get_parameters(thapca_flat, iterations = 10^5)

# Summarizes the flat prior analysis
thapca_flat_summary <- summary(thapca_flat)
mean_eff <- thapca_flat_summary$fixed[2, 1]
mean_sd <- thapca_flat_summary$fixed[2, 2]

# Some summaries used for plot legend
any_harm <- sum(posteriors$b_TreatRand > 0) / nrow(posteriors)
severe_harm <- sum(posteriors$b_TreatRand > log(1.25)) / nrow(posteriors)
rope <- sum(posteriors$b_TreatRand < log(1.1) &
            posteriors$b_TreatRand > log(1/1.1)) /
        nrow(posteriors)

# Figure 1
temp_plot <- ggplot(posteriors, aes(x = b_TreatRand)) +
             geom_density()
p <- ggplot_build(temp_plot)
plot_data <- p$data[[1]][, c(1, 2)]

fig1 <- ggplot(plot_data, aes(x = x, y = y, group = 1)) +
  geom_line(size = 1) +
  geom_area(mapping = aes(x = ifelse(x > 0, x, 0)),
            fill = "#FFB266") +
  geom_area(mapping = aes(x = ifelse(x < 0, x, 0)),
            fill = "lightblue") +
  geom_area(mapping = aes(x = ifelse(x > log(1.25),
                                     x, 0)),
            fill = "#FF9933") +
  geom_segment(inherit.aes = FALSE,
               data =
                 subset(plot_data,
                        x > log(1 / 1.1) &
                        x < log(1.1))[seq(1, nrow(subset(plot_data,
                                                         x > log(1/1.1) &
                                                          x < log(1.1))), 5), ],
               aes(x = x, y = 0, xend = x, yend = y)) +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1) +
  geom_vline(xintercept = log(1.25),
             color = "#FF9933", linetype = 1) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_classic() +
  labs(x = "log(OR)", y = "") +
  scale_x_continuous(expand = c(0, 0),
                     labels = seq(-2, 0.6, by = 0.2),
                     breaks = seq(-2, 0.6, by = 0.2)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(-2.1, 0.7),
                  ylim = c(0, 1.5)) +
  annotate("text", x = 0.35, y = 1.3,
           label = paste0("P(severe harm) = ",
                          round(severe_harm, 2), " ",
                          sprintf('\u21d2')))+
  annotate("text", x = 0.1, y = 1.3,
           label = paste0("P(harm) = ",
                          round(any_harm, 2), " ",
                          sprintf('\u21d2')))+
  annotate("text", x = -1, y = 1.3,
           label = paste0(sprintf('\u21d0'), " ",
                          "P(benefit) = ",
                          1 - round(any_harm, 2))) +
  #  annotate("rect", xmin = -0.06, xmax = 0.06, ymin = 0.06,
             # ymax = 0.18,fill="white")+
  annotate("label", x = 0, y = 0.1,
           label =
             paste0("ROPE =",
                    round(sum(posteriors$b_TreatRand < log(1.1) &
                                posteriors$b_TreatRand > log(1/1.1)) /
                            nrow(posteriors), 2)))

# Define priors
# Helper functions to determine priors
getDiff <- function(x, priormean = 0,
                    propbelow, belowcutoff) {
  propbelow - pnorm(belowcutoff, priormean, x)
}
getPriorSD <- function(priormean = 0,
                       propbelow, belowcutoff) {
  
  root <- uniroot(getDiff, propbelow = propbelow,
                  belowcutoff = belowcutoff,
                  priormean = priormean,
                  interval = c(0, 10))$root
  return(round(root*200) / 200)
  
}

# Neutral moderate prior (SD = 0.355)
getPriorSD(propbelow = 0.025, belowcutoff = log(0.5))
m1priors <- prior(normal(0, 0.355), class = "b")
m1 <- brm(PrimaryEndpoint ~ TreatRand + AgeGroup,
          data = dat_primary, prior = m1priors,
          family = "bernoulli", seed = 123)
posteriors_m1 <- insight::get_parameters(m1)[, 2]
any_harm_m1 <- sum(posteriors_m1 > 0) / length(posteriors_m1)
severe_harm_m1 <- sum(posteriors_m1 > log(1.25)) / length(posteriors_m1)

# Optimistic moderate prior (SD = 0.400)
getPriorSD(propbelow = 0.85, belowcutoff = 0, priormean = log(0.66))
m2priors <- prior(normal(log(0.66), 0.4), class = "b")
m2 <- brm(PrimaryEndpoint ~ TreatRand + AgeGroup,
          data = dat_primary, prior = m2priors,
          family = "bernoulli", seed = 123)
posteriors_m2 <- insight::get_parameters(m2)[, 2]
any_harm_m2 <- sum(posteriors_m2 > 0) / length(posteriors_m2)
severe_harm_m2 <- sum(posteriors_m2 > log(1.25)) / length(posteriors_m2)

# Pessimistic weak prior (SD = 0.790)
getPriorSD(propbelow = 0.3, belowcutoff = 0, priormean = -log(0.66))
m3priors <- prior(normal(-log(0.66), 0.79), class = "b")
m3 <- brm(PrimaryEndpoint ~ TreatRand + AgeGroup,
          data = dat_primary, prior = m3priors,
          family = "bernoulli", seed = 123)
posteriors_m3 <- insight::get_parameters(m3)[, 2]
any_harm_m3 <- sum(posteriors_m3 > 0) / length(posteriors_m3)
severe_harm_m3 <- sum(posteriors_m3 > log(1.25)) / length(posteriors_m3)

#######################################################################
# Bayesian re-analysis of secondary outcome 1

# Flat prior
thapca_flat2 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
                    data = dat_secondary1,
                    prior = flat_prior,
                    family = "bernoulli", seed = 123)
posteriors <- insight::get_parameters(thapca_flat2, iterations = 10^5)

# Summarizes the flat prior analysis
summary(thapca_flat2)

# Some summaries used for plot legend
any_harm <- sum(posteriors$b_TreatRand > 0) / nrow(posteriors)
severe_harm <- sum(posteriors$b_TreatRand > log(1.25)) / nrow(posteriors)
rope <- sum(posteriors$b_TreatRand < log(1.1) &
              posteriors$b_TreatRand > log(1/1.1)) /
  nrow(posteriors)

# Neutral moderate prior (SD = 0.355)
m4 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m1priors,
          family = "bernoulli", seed = 123)
posteriors_m4 <- insight::get_parameters(m4)[, 2]
any_harm_m4 <- sum(posteriors_m4 > 0) / length(posteriors_m4)
severe_harm_m4 <- sum(posteriors_m4 > log(1.25)) / length(posteriors_m4)

# Optimistic moderate prior (SD = 0.400)
m5 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m2priors,
          family = "bernoulli", seed = 123)
posteriors_m5 <- insight::get_parameters(m5)[, 2]
any_harm_m5 <- sum(posteriors_m5 > 0) / length(posteriors_m5)
severe_harm_m5 <- sum(posteriors_m5 > log(1.25)) / length(posteriors_m5)

# Pessimistic weak prior (SD = 0.790)
m6 <- brm(SurviveM12 ~ TreatRand + AgeGroup,
          data = dat_secondary1, prior = m3priors,
          family = "bernoulli", seed = 123)
posteriors_m6 <- insight::get_parameters(m6)[, 2]
any_harm_m6 <- sum(posteriors_m6 > 0) / length(posteriors_m6)
severe_harm_m6 <- sum(posteriors_m6 > log(1.25)) / length(posteriors_m6)
