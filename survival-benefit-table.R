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
library(rstanarm)

# Load trial data
dat <- read.csv("../outcomes.csv")
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)

# Frequentist RD model
mod_freq <- glm(PrimaryEndpoint ~ Trt + AgeGroup,
                data = dat_primary,
                family = binomial(link = "identity"))

mod_freq2 <- glm(PrimaryEndpoint ~ Trt + AgeGroup,
                 data = dat_primary,
                 family = gaussian(link = "identity"))

# Bayesian analysis
# One prior for now (Moderate optimistic)
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = 0.1)
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 5),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(mo_sd, 5)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

mean(posteriors_mo$Trt > 0)
mean(posteriors_mo$Trt > 0.02)
mean(posteriors_mo$Trt > 0.04)
mean(posteriors_mo$Trt > 0.06)
mean(posteriors_mo$Trt > 0.08)
mean(posteriors_mo$Trt > 0.1)
mean(posteriors_mo$Trt > 0.2)
