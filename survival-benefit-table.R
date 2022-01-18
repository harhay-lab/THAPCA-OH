# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(flextable)
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

#######################################################################
# First do analysis on risk difference scale

# Optimistic priors
# Strong optimistic prior
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = 0.1)
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(so_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)

so_med <- median(posteriors_so$Trt)
so_low <- quantile(posteriors_so$Trt, 0.025)
so_up <- quantile(posteriors_so$Trt, 0.975)

so0.00 <- mean(posteriors_so$Trt > 0)
so0.02 <- mean(posteriors_so$Trt > 0.02)
so0.05 <- mean(posteriors_so$Trt > 0.05)
so0.10 <- mean(posteriors_so$Trt > 0.1)
so0.15 <- mean(posteriors_so$Trt > 0.15)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = 0.1)
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(mo_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

mo_med <- median(posteriors_mo$Trt)
mo_low <- quantile(posteriors_mo$Trt, 0.025)
mo_up <- quantile(posteriors_mo$Trt, 0.975)

mo0.00 <- mean(posteriors_mo$Trt > 0)
mo0.02 <- mean(posteriors_mo$Trt > 0.02)
mo0.05 <- mean(posteriors_mo$Trt > 0.05)
mo0.10 <- mean(posteriors_mo$Trt > 0.1)
mo0.15 <- mean(posteriors_mo$Trt > 0.15)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = 0.1)
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0.1, 0),
                                  scale = c(wo_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wo <- insight::get_parameters(wo_mod)

wo_med <- median(posteriors_wo$Trt)
wo_low <- quantile(posteriors_wo$Trt, 0.025)
wo_up <- quantile(posteriors_wo$Trt, 0.975)

wo0.00 <- mean(posteriors_wo$Trt > 0)
wo0.02 <- mean(posteriors_wo$Trt > 0.02)
wo0.05 <- mean(posteriors_wo$Trt > 0.05)
wo0.10 <- mean(posteriors_wo$Trt > 0.1)
wo0.15 <- mean(posteriors_wo$Trt > 0.15)

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = -0.05,
                    priormean = 0)
sn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(sn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_sn <- insight::get_parameters(sn_mod)

sn_med <- median(posteriors_sn$Trt)
sn_low <- quantile(posteriors_sn$Trt, 0.025)
sn_up <- quantile(posteriors_sn$Trt, 0.975)

sn0.00 <- mean(posteriors_sn$Trt > 0)
sn0.02 <- mean(posteriors_sn$Trt > 0.02)
sn0.05 <- mean(posteriors_sn$Trt > 0.05)
sn0.10 <- mean(posteriors_sn$Trt > 0.1)
sn0.15 <- mean(posteriors_sn$Trt > 0.15)

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = -0.1,
                    priormean = 0)
mn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mn <- insight::get_parameters(mn_mod)

mn_med <- median(posteriors_mn$Trt)
mn_low <- quantile(posteriors_mn$Trt, 0.025)
mn_up <- quantile(posteriors_mn$Trt, 0.975)

mn0.00 <- mean(posteriors_mn$Trt > 0)
mn0.02 <- mean(posteriors_mn$Trt > 0.02)
mn0.05 <- mean(posteriors_mn$Trt > 0.05)
mn0.10 <- mean(posteriors_mn$Trt > 0.1)
mn0.15 <- mean(posteriors_mn$Trt > 0.15)

# Weak neutral prior
wn_sd <- 5
wn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(wn_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wn <- insight::get_parameters(wn_mod)

wn_med <- median(posteriors_wn$Trt)
wn_low <- quantile(posteriors_wn$Trt, 0.025)
wn_up <- quantile(posteriors_wn$Trt, 0.975)

wn0.00 <- mean(posteriors_wn$Trt > 0)
wn0.02 <- mean(posteriors_wn$Trt > 0.02)
wn0.05 <- mean(posteriors_wn$Trt > 0.05)
wn0.10 <- mean(posteriors_wn$Trt > 0.1)
wn0.15 <- mean(posteriors_wn$Trt > 0.15)

# Pessimistic priors
# Strong pessimistic prior
sp_sd <- so_sd
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(-0.1, 0),
                                  scale = c(sp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)

sp_med <- median(posteriors_sp$Trt)
sp_low <- quantile(posteriors_sp$Trt, 0.025)
sp_up <- quantile(posteriors_sp$Trt, 0.975)

sp0.00 <- mean(posteriors_sp$Trt > 0)
sp0.02 <- mean(posteriors_sp$Trt > 0.02)
sp0.05 <- mean(posteriors_sp$Trt > 0.05)
sp0.10 <- mean(posteriors_sp$Trt > 0.1)
sp0.15 <- mean(posteriors_sp$Trt > 0.15)

# Moderate pessimistic prior
mp_sd <- mo_sd
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)

mp_med <- median(posteriors_mp$Trt)
mp_low <- quantile(posteriors_mp$Trt, 0.025)
mp_up <- quantile(posteriors_mp$Trt, 0.975)

mp0.00 <- mean(posteriors_mp$Trt > 0)
mp0.02 <- mean(posteriors_mp$Trt > 0.02)
mp0.05 <- mean(posteriors_mp$Trt > 0.05)
mp0.10 <- mean(posteriors_mp$Trt > 0.1)
mp0.15 <- mean(posteriors_mp$Trt > 0.15)

# Weak pessimistic prior
wp_sd <- wo_sd
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(wp_sd, 100)),
                   family = gaussian(link = "identity"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)

wp_med <- median(posteriors_wp$Trt)
wp_low <- quantile(posteriors_wp$Trt, 0.025)
wp_up <- quantile(posteriors_wp$Trt, 0.975)

wp0.00 <- mean(posteriors_wp$Trt > 0)
wp0.02 <- mean(posteriors_wp$Trt > 0.02)
wp0.05 <- mean(posteriors_wp$Trt > 0.05)
wp0.10 <- mean(posteriors_wp$Trt > 0.1)
wp0.15 <- mean(posteriors_wp$Trt > 0.15)


# Make table
df <- data.frame("Prior Belief" = c("Optimistic", "", "",
                                    "Neutral", "", "",
                                    "Pessimistic", "", ""),
                 "Strength" = rep(c("Weak", "Moderate", "Strong"), 3),
                 "Median" = c(paste(round(wo_med, 2), " (",
                                    round(wo_low, 2), ", ",
                                    round(wo_up, 2), ")", sep = ""),
                              paste(round(mo_med, 2), " (",
                                    round(mo_low, 2), ", ",
                                    round(mo_up, 2), ")", sep = ""),
                              paste(round(so_med, 2), " (",
                                    round(so_low, 2), ", ",
                                    round(so_up, 2), ")", sep = ""),
                              paste(round(wn_med, 2), " (",
                                    round(wn_low, 2), ", ",
                                    round(wn_up, 2), ")", sep = ""),
                              paste(round(mn_med, 2), " (",
                                    round(mn_low, 2), ", ",
                                    round(mn_up, 2), ")", sep = ""),
                              paste(round(sn_med, 2), " (",
                                    round(sn_low, 2), ", ",
                                    round(sn_up, 2), ")", sep = ""),
                              paste(round(wp_med, 2), " (",
                                    round(wp_low, 2), ", ",
                                    round(wp_up, 2), ")", sep = ""),
                              paste(round(mp_med, 2), " (",
                                    round(mp_low, 2), ", ",
                                    round(mp_up, 2), ")", sep = ""),
                              paste(round(sp_med, 2), " (",
                                    round(sp_low, 2), ", ",
                                    round(sp_up, 2), ")", sep = "")),
                 ">0" = 100*round(c(wo0.00, mo0.00, so0.00,
                                    wn0.00, mn0.00, sn0.00,
                                    wp0.00, mp0.00, sp0.00), 2),
                 ">0.02" = 100*round(c(wo0.02, mo0.02, so0.02,
                                       wn0.02, mn0.02, sn0.02,
                                       wp0.02, mp0.02, sp0.02), 2),
                 ">0.05" = 100*round(c(wo0.05, mo0.05, so0.05,
                                       wn0.05, mn0.05, sn0.05,
                                       wp0.05, mp0.05, sp0.05), 2),
                 ">0.1" = 100*round(c(wo0.10, mo0.10, so0.10,
                                      wn0.10, mn0.10, sn0.10,
                                      wp0.10, mp0.10, sp0.10), 2),
                 ">0.15" = 100*round(c(wo0.15, mo0.15, so0.15,
                                       wn0.15, mn0.15, sn0.15,
                                       wp0.15, mp0.15, sp0.15), 2))

ft <- flextable(df) %>%
  theme_zebra()
ft <- width(ft, j = 3, width = 1.5)

# Export table to Word
read_docx() %>% 
  body_add_par("Table 2 (RD Scale)") %>% 
  body_add_flextable(value = ft) %>% 
  print(target = "example_table_word.docx")


#######################################################################
# Now do analysis on relative risk scale

# Optimistic priors
# Strong optimistic prior
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = log(1.25))
so_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1.25), 0),
                                  scale = c(so_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_so <- insight::get_parameters(so_mod)

so_med <- median(exp(posteriors_so$Trt))
so_low <- quantile(exp(posteriors_so$Trt), 0.025)
so_up <- quantile(exp(posteriors_so$Trt), 0.975)

so1.00 <- mean(exp(posteriors_so$Trt) > 1)
so1.10 <- mean(exp(posteriors_so$Trt) > 1.1)
so1.25 <- mean(exp(posteriors_so$Trt) > 1.25)
so1.50 <- mean(exp(posteriors_so$Trt) > 1.5)
so2.00 <- mean(exp(posteriors_so$Trt) > 2)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = log(1.25))
mo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1.25), 0),
                                  scale = c(mo_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mo <- insight::get_parameters(mo_mod)

mo_med <- median(exp(posteriors_mo$Trt))
mo_low <- quantile(exp(posteriors_mo$Trt), 0.025)
mo_up <- quantile(exp(posteriors_mo$Trt), 0.975)

mo1.00 <- mean(exp(posteriors_mo$Trt) > 1)
mo1.10 <- mean(exp(posteriors_mo$Trt) > 1.1)
mo1.25 <- mean(exp(posteriors_mo$Trt) > 1.25)
mo1.50 <- mean(exp(posteriors_mo$Trt) > 1.5)
mo2.00 <- mean(exp(posteriors_mo$Trt) > 2)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = log(1.25))
wo_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1.25), 0),
                                  scale = c(wo_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wo <- insight::get_parameters(wo_mod)

wo_med <- median(exp(posteriors_wo$Trt))
wo_low <- quantile(exp(posteriors_wo$Trt), 0.025)
wo_up <- quantile(exp(posteriors_wo$Trt), 0.975)

wo1.00 <- mean(exp(posteriors_wo$Trt) > 1)
wo1.10 <- mean(exp(posteriors_wo$Trt) > 1.1)
wo1.25 <- mean(exp(posteriors_wo$Trt) > 1.25)
wo1.50 <- mean(exp(posteriors_wo$Trt) > 1.5)
wo2.00 <- mean(exp(posteriors_wo$Trt) > 2)

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5))
sn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(sn_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sn <- insight::get_parameters(sn_mod)

sn_med <- median(exp(posteriors_sn$Trt))
sn_low <- quantile(exp(posteriors_sn$Trt), 0.025)
sn_up <- quantile(exp(posteriors_sn$Trt), 0.975)

sn1.00 <- mean(exp(posteriors_sn$Trt) > 1)
sn1.10 <- mean(exp(posteriors_sn$Trt) > 1.1)
sn1.25 <- mean(exp(posteriors_sn$Trt) > 1.25)
sn1.50 <- mean(exp(posteriors_sn$Trt) > 1.5)
sn2.00 <- mean(exp(posteriors_sn$Trt) > 2)

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5))
mn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(mn_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mn <- insight::get_parameters(mn_mod)

mn_med <- median(exp(posteriors_mn$Trt))
mn_low <- quantile(exp(posteriors_mn$Trt), 0.025)
mn_up <- quantile(exp(posteriors_mn$Trt), 0.975)

mn1.00 <- mean(exp(posteriors_mn$Trt) > 1)
mn1.10 <- mean(exp(posteriors_mn$Trt) > 1.1)
mn1.25 <- mean(exp(posteriors_mn$Trt) > 1.25)
mn1.50 <- mean(exp(posteriors_mn$Trt) > 1.5)
mn2.00 <- mean(exp(posteriors_mn$Trt) > 2)

# Weak neutral prior (SD = 5)
wn_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(0, 0),
                                  scale = c(5, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wn <- insight::get_parameters(wn_mod)

wn_med <- median(exp(posteriors_wn$Trt))
wn_low <- quantile(exp(posteriors_wn$Trt), 0.025)
wn_up <- quantile(exp(posteriors_wn$Trt), 0.975)

wn1.00 <- mean(exp(posteriors_wn$Trt) > 1)
wn1.10 <- mean(exp(posteriors_wn$Trt) > 1.1)
wn1.25 <- mean(exp(posteriors_wn$Trt) > 1.25)
wn1.50 <- mean(exp(posteriors_wn$Trt) > 1.5)
wn2.00 <- mean(exp(posteriors_wn$Trt) > 2)

# Pessimistic priors
# Strong pessimistic prior
sp_sd <- so_sd
sp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1/1.25), 0),
                                  scale = c(sp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_sp <- insight::get_parameters(sp_mod)

sp_med <- median(exp(posteriors_sp$Trt))
sp_low <- quantile(exp(posteriors_sp$Trt), 0.025)
sp_up <- quantile(exp(posteriors_sp$Trt), 0.975)

sp1.00 <- mean(exp(posteriors_sp$Trt) > 1)
sp1.10 <- mean(exp(posteriors_sp$Trt) > 1.1)
sp1.25 <- mean(exp(posteriors_sp$Trt) > 1.25)
sp1.50 <- mean(exp(posteriors_sp$Trt) > 1.5)
sp2.00 <- mean(exp(posteriors_sp$Trt) > 2)

# Moderate pessimistic prior
mp_sd <- mo_sd
mp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1/1.25), 0),
                                  scale = c(mp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_mp <- insight::get_parameters(mp_mod)

mp_med <- median(exp(posteriors_mp$Trt))
mp_low <- quantile(exp(posteriors_mp$Trt), 0.025)
mp_up <- quantile(exp(posteriors_mp$Trt), 0.975)

mp1.00 <- mean(exp(posteriors_mp$Trt) > 1)
mp1.10 <- mean(exp(posteriors_mp$Trt) > 1.1)
mp1.25 <- mean(exp(posteriors_mp$Trt) > 1.25)
mp1.50 <- mean(exp(posteriors_mp$Trt) > 1.5)
mp2.00 <- mean(exp(posteriors_mp$Trt) > 2)

# Weak pessimistic prior
wp_sd <- wo_sd
wp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeGroup,
                   data = dat_primary,
                   prior_intercept = normal(0, 100),
                   prior = normal(location = c(log(1/1.25), 0),
                                  scale = c(wp_sd, 100)),
                   family = binomial(link = "log"), seed = 1234)
posteriors_wp <- insight::get_parameters(wp_mod)

wp_med <- median(exp(posteriors_wp$Trt))
wp_low <- quantile(exp(posteriors_wp$Trt), 0.025)
wp_up <- quantile(exp(posteriors_wp$Trt), 0.975)

wp1.00 <- mean(exp(posteriors_wp$Trt) > 1)
wp1.10 <- mean(exp(posteriors_wp$Trt) > 1.1)
wp1.25 <- mean(exp(posteriors_wp$Trt) > 1.25)
wp1.50 <- mean(exp(posteriors_wp$Trt) > 1.5)
wp2.00 <- mean(exp(posteriors_wp$Trt) > 2)


# Make table
df2 <- data.frame("Prior Belief" = c("Optimistic", "", "",
                                    "Neutral", "", "",
                                    "Pessimistic", "", ""),
                 "Strength" = rep(c("Weak", "Moderate", "Strong"), 3),
                 "Median" = c(paste(round(wo_med, 2), " (",
                                    round(wo_low, 2), ", ",
                                    round(wo_up, 2), ")", sep = ""),
                              paste(round(mo_med, 2), " (",
                                    round(mo_low, 2), ", ",
                                    round(mo_up, 2), ")", sep = ""),
                              paste(round(so_med, 2), " (",
                                    round(so_low, 2), ", ",
                                    round(so_up, 2), ")", sep = ""),
                              paste(round(wn_med, 2), " (",
                                    round(wn_low, 2), ", ",
                                    round(wn_up, 2), ")", sep = ""),
                              paste(round(mn_med, 2), " (",
                                    round(mn_low, 2), ", ",
                                    round(mn_up, 2), ")", sep = ""),
                              paste(round(sn_med, 2), " (",
                                    round(sn_low, 2), ", ",
                                    round(sn_up, 2), ")", sep = ""),
                              paste(round(wp_med, 2), " (",
                                    round(wp_low, 2), ", ",
                                    round(wp_up, 2), ")", sep = ""),
                              paste(round(mp_med, 2), " (",
                                    round(mp_low, 2), ", ",
                                    round(mp_up, 2), ")", sep = ""),
                              paste(round(sp_med, 2), " (",
                                    round(sp_low, 2), ", ",
                                    round(sp_up, 2), ")", sep = "")),
                 ">1" = 100*round(c(wo1.00, mo1.00, so1.00,
                                    wn1.00, mn1.00, sn1.00,
                                    wp1.00, mp1.00, sp1.00), 2),
                 ">1.1" = 100*round(c(wo1.10, mo1.10, so1.10,
                                      wn1.10, mn1.10, sn1.10,
                                      wp1.10, mp1.10, sp1.10), 2),
                 ">1.25" = 100*round(c(wo1.25, mo1.25, so1.25,
                                       wn1.25, mn1.25, sn1.25,
                                       wp1.25, mp1.25, sp1.25), 2),
                 ">1.5" = 100*round(c(wo1.50, mo1.50, so1.50,
                                      wn1.50, mn1.50, sn1.50,
                                      wp1.50, mp1.50, sp1.50), 2),
                 ">2" = 100*round(c(wo2.00, mo2.00, so2.00,
                                    wn2.00, mn2.00, sn2.00,
                                    wp2.00, mp2.00, sp2.00), 2))

ft2 <- flextable(df2) %>%
  theme_zebra()
ft2 <- width(ft2, j = 3, width = 1.5)

# Export table to Word
read_docx() %>% 
  body_add_par("Table 2 (RD Scale)") %>% 
  body_add_flextable(value = ft) %>% 
  print(target = "example_table_word.docx")