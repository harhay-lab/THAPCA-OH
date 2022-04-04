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
library(multipanelfigure)
library(pvaluefunctions)
library(sanon)
#library(rstanarm)

# Load trial data
dat <- read.csv("/Users/blette/Downloads/outcomes.csv")
dat_secondary1 <- dat[dat$SurviveM12 != 97, ]
dat_secondary1$Trt <- abs(2 - dat_secondary1$TreatRand)
dat_secondary1$AgeFactor <- as.factor(dat_secondary1$AgeGroup)

# Helper functions to determine priors for Bayesian analysis
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
# First run analysis and get RD and RR results

# Optimistic priors
# Strong optimistic prior
o_mean <- log(1.25)
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(so_sd)
so_prior <- prior(normal(o_mean, so_sd), class = "b", coef = Trt)
so_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = so_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_so <- posterior_epred(so_mod,
                            newdata = mutate(so_mod$data, Trt = 1))
pred0_so <- posterior_epred(so_mod,
                            newdata = mutate(so_mod$data, Trt = 0))

ratio_so <- apply(pred1_so, 1, mean) / apply(pred0_so, 1, mean)
diff_so <- 100*(apply(pred1_so, 1, mean) - apply(pred0_so, 1, mean))

so_med_d <- median(diff_so)
so_low_d <- quantile(diff_so, 0.025)
so_up_d <- quantile(diff_so, 0.975)
so0.00 <- mean(diff_so > 0)
so0.02 <- mean(diff_so > 2)
so0.05 <- mean(diff_so > 5)
so0.10 <- mean(diff_so > 10)
so0.15 <- mean(diff_so > 15)

so_med_r <- median(ratio_so)
so_low_r <- quantile(ratio_so, 0.025)
so_up_r <- quantile(ratio_so, 0.975)
so1.00 <- mean(ratio_so > 1)
so1.10 <- mean(ratio_so > 1.1)
so1.25 <- mean(ratio_so > 1.25)
so1.50 <- mean(ratio_so > 1.5)
so2.00 <- mean(ratio_so > 2)

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(mo_sd)
mo_prior <- prior(normal(o_mean, mo_sd), class = "b", coef = Trt)
mo_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = mo_prior, stanvars = stanvars, iter = 10000,
              family = "bernoulli", seed = 1234, chains = 8)

pred1_mo <- posterior_epred(mo_mod,
                            newdata = mutate(mo_mod$data,Trt = 1))
pred0_mo <- posterior_epred(mo_mod,
                            newdata = mutate(mo_mod$data, Trt = 0))

ratio_mo <- apply(pred1_mo, 1, mean) / apply(pred0_mo, 1, mean)
diff_mo <- 100*(apply(pred1_mo, 1, mean) - apply(pred0_mo, 1, mean))

mo_med_d <- median(diff_mo)
mo_low_d <- quantile(diff_mo, 0.025)
mo_up_d <- quantile(diff_mo, 0.975)
mo0.00 <- mean(diff_mo > 0)
mo0.02 <- mean(diff_mo > 2)
mo0.05 <- mean(diff_mo > 5)
mo0.10 <- mean(diff_mo > 10)
mo0.15 <- mean(diff_mo > 15)

mo_med_r <- median(ratio_mo)
mo_low_r <- quantile(ratio_mo, 0.025)
mo_up_r <- quantile(ratio_mo, 0.975)
mo1.00 <- mean(ratio_mo > 1)
mo1.10 <- mean(ratio_mo > 1.1)
mo1.25 <- mean(ratio_mo > 1.25)
mo1.50 <- mean(ratio_mo > 1.5)
mo2.00 <- mean(ratio_mo > 2)

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(wo_sd)
wo_prior <- prior(normal(o_mean, wo_sd), class = "b", coef = Trt)
wo_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = wo_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wo <- posterior_epred(wo_mod,
                            newdata = mutate(wo_mod$data, Trt = 1))
pred0_wo <- posterior_epred(wo_mod,
                            newdata = mutate(wo_mod$data, Trt = 0))

ratio_wo <- apply(pred1_wo, 1, mean) / apply(pred0_wo, 1, mean)
diff_wo <- 100*(apply(pred1_wo, 1, mean) - apply(pred0_wo, 1, mean))

wo_med_d <- median(diff_wo)
wo_low_d <- quantile(diff_wo, 0.025)
wo_up_d <- quantile(diff_wo, 0.975)
wo0.00 <- mean(diff_wo > 0)
wo0.02 <- mean(diff_wo > 2)
wo0.05 <- mean(diff_wo > 5)
wo0.10 <- mean(diff_wo > 10)
wo0.15 <- mean(diff_wo > 15)

wo_med_r <- median(ratio_wo)
wo_low_r <- quantile(ratio_wo, 0.025)
wo_up_r <- quantile(ratio_wo, 0.975)
wo1.00 <- mean(ratio_wo > 1)
wo1.10 <- mean(ratio_wo > 1.1)
wo1.25 <- mean(ratio_wo > 1.25)
wo1.50 <- mean(ratio_wo > 1.5)
wo2.00 <- mean(ratio_wo > 2)

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5),
                    priormean = 0)
stanvars <- stanvar(sn_sd)
sn_prior <- prior(normal(0, sn_sd), class = "b", coef = Trt)
sn_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = sn_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_sn <- posterior_epred(sn_mod,
                            newdata = mutate(sn_mod$data, Trt = 1))
pred0_sn <- posterior_epred(sn_mod,
                            newdata = mutate(sn_mod$data, Trt = 0))

ratio_sn <- apply(pred1_sn, 1, mean) / apply(pred0_sn, 1, mean)
diff_sn <- 100*(apply(pred1_sn, 1, mean) - apply(pred0_sn, 1, mean))

sn_med_d <- median(diff_sn)
sn_low_d <- quantile(diff_sn, 0.025)
sn_up_d <- quantile(diff_sn, 0.975)
sn0.00 <- mean(diff_sn > 0)
sn0.02 <- mean(diff_sn > 2)
sn0.05 <- mean(diff_sn > 5)
sn0.10 <- mean(diff_sn > 10)
sn0.15 <- mean(diff_sn > 15)

sn_med_r <- median(ratio_sn)
sn_low_r <- quantile(ratio_sn, 0.025)
sn_up_r <- quantile(ratio_sn, 0.975)
sn1.00 <- mean(ratio_sn > 1)
sn1.10 <- mean(ratio_sn > 1.1)
sn1.25 <- mean(ratio_sn > 1.25)
sn1.50 <- mean(ratio_sn > 1.5)
sn2.00 <- mean(ratio_sn > 2)

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5),
                    priormean = 0)
stanvars <- stanvar(mn_sd)
mn_prior <- prior(normal(0, mn_sd), class = "b", coef = Trt)
mn_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = mn_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_mn <- posterior_epred(mn_mod,
                            newdata = mutate(mn_mod$data, Trt = 1))
pred0_mn <- posterior_epred(mn_mod,
                            newdata = mutate(mn_mod$data, Trt = 0))

ratio_mn <- apply(pred1_mn, 1, mean) / apply(pred0_mn, 1, mean)
diff_mn <- 100*(apply(pred1_mn, 1, mean) - apply(pred0_mn, 1, mean))

mn_med_d <- median(diff_mn)
mn_low_d <- quantile(diff_mn, 0.025)
mn_up_d <- quantile(diff_mn, 0.975)
mn0.00 <- mean(diff_mn > 0)
mn0.02 <- mean(diff_mn > 2)
mn0.05 <- mean(diff_mn > 5)
mn0.10 <- mean(diff_mn > 10)
mn0.15 <- mean(diff_mn > 15)

mn_med_r <- median(ratio_mn)
mn_low_r <- quantile(ratio_mn, 0.025)
mn_up_r <- quantile(ratio_mn, 0.975)
mn1.00 <- mean(ratio_mn > 1)
mn1.10 <- mean(ratio_mn > 1.1)
mn1.25 <- mean(ratio_mn > 1.25)
mn1.50 <- mean(ratio_mn > 1.5)
mn2.00 <- mean(ratio_mn > 2)

# Weak neutral prior
wn_sd <- 3
wn_prior <- prior(normal(0, 3), class = "b", coef = Trt)
wn_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = wn_prior, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wn <- posterior_epred(wn_mod,
                            newdata = mutate(wn_mod$data, Trt = 1))
pred0_wn <- posterior_epred(wn_mod,
                            newdata = mutate(wn_mod$data, Trt = 0))

ratio_wn <- apply(pred1_wn, 1, mean) / apply(pred0_wn, 1, mean)
diff_wn <- 100*(apply(pred1_wn, 1, mean) - apply(pred0_wn, 1, mean))

wn_med_d <- median(diff_wn)
wn_low_d <- quantile(diff_wn, 0.025)
wn_up_d <- quantile(diff_wn, 0.975)
wn0.00 <- mean(diff_wn > 0)
wn0.02 <- mean(diff_wn > 2)
wn0.05 <- mean(diff_wn > 5)
wn0.10 <- mean(diff_wn > 10)
wn0.15 <- mean(diff_wn > 15)

wn_med_r <- median(ratio_wn)
wn_low_r <- quantile(ratio_wn, 0.025)
wn_up_r <- quantile(ratio_wn, 0.975)
wn1.00 <- mean(ratio_wn > 1)
wn1.10 <- mean(ratio_wn > 1.1)
wn1.25 <- mean(ratio_wn > 1.25)
wn1.50 <- mean(ratio_wn > 1.5)
wn2.00 <- mean(ratio_wn > 2)

# Pessimistic priors
# Strong pessimistic prior
p_mean <- log(1/1.25)
sp_sd <- so_sd
stanvars <- stanvar(p_mean) + stanvar(sp_sd)
sp_prior <- prior(normal(p_mean, sp_sd), class = "b", coef = Trt)
sp_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = sp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_sp <- posterior_epred(sp_mod,
                            newdata = mutate(sp_mod$data, Trt = 1))
pred0_sp <- posterior_epred(sp_mod,
                            newdata = mutate(sp_mod$data, Trt = 0))

ratio_sp <- apply(pred1_sp, 1, mean) / apply(pred0_sp, 1, mean)
diff_sp <- 100*(apply(pred1_sp, 1, mean) - apply(pred0_sp, 1, mean))

sp_med_d <- median(diff_sp)
sp_low_d <- quantile(diff_sp, 0.025)
sp_up_d <- quantile(diff_sp, 0.975)
sp0.00 <- mean(diff_sp > 0)
sp0.02 <- mean(diff_sp > 2)
sp0.05 <- mean(diff_sp > 5)
sp0.10 <- mean(diff_sp > 10)
sp0.15 <- mean(diff_sp > 15)

sp_med_r <- median(ratio_sp)
sp_low_r <- quantile(ratio_sp, 0.025)
sp_up_r <- quantile(ratio_sp, 0.975)
sp1.00 <- mean(ratio_sp > 1)
sp1.10 <- mean(ratio_sp > 1.1)
sp1.25 <- mean(ratio_sp > 1.25)
sp1.50 <- mean(ratio_sp > 1.5)
sp2.00 <- mean(ratio_sp > 2)

# Moderate pessimistic prior
mp_sd <- mo_sd
stanvars <- stanvar(p_mean) + stanvar(mp_sd)
mp_prior <- prior(normal(p_mean, mp_sd), class = "b", coef = Trt)
mp_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = mp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_mp <- posterior_epred(mp_mod,
                            newdata = mutate(mp_mod$data, Trt = 1))
pred0_mp <- posterior_epred(mp_mod,
                            newdata = mutate(mp_mod$data, Trt = 0))

ratio_mp <- apply(pred1_mp, 1, mean) / apply(pred0_mp, 1, mean)
diff_mp <- 100*(apply(pred1_mp, 1, mean) - apply(pred0_mp, 1, mean))

mp_med_d <- median(diff_mp)
mp_low_d <- quantile(diff_mp, 0.025)
mp_up_d <- quantile(diff_mp, 0.975)
mp0.00 <- mean(diff_mp > 0)
mp0.02 <- mean(diff_mp > 2)
mp0.05 <- mean(diff_mp > 5)
mp0.10 <- mean(diff_mp > 10)
mp0.15 <- mean(diff_mp > 15)

mp_med_r <- median(ratio_mp)
mp_low_r <- quantile(ratio_mp, 0.025)
mp_up_r <- quantile(ratio_mp, 0.975)
mp1.00 <- mean(ratio_mp > 1)
mp1.10 <- mean(ratio_mp > 1.1)
mp1.25 <- mean(ratio_mp > 1.25)
mp1.50 <- mean(ratio_mp > 1.5)
mp2.00 <- mean(ratio_mp > 2)

# Weak pessimistic prior
wp_sd <- wo_sd
stanvars <- stanvar(p_mean) + stanvar(wp_sd)
wp_prior <- prior(normal(p_mean, wp_sd), class = "b", coef = Trt)
wp_mod <- brm(SurviveM12 ~ Trt + AgeFactor, data = dat_secondary1,
              prior = wp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wp <- posterior_epred(wp_mod,
                            newdata = mutate(wp_mod$data, Trt = 1))
pred0_wp <- posterior_epred(wp_mod,
                            newdata = mutate(wp_mod$data, Trt = 0))

ratio_wp <- apply(pred1_wp, 1, mean) / apply(pred0_wp, 1, mean)
diff_wp <- 100*(apply(pred1_wp, 1, mean) - apply(pred0_wp, 1, mean))

wp_med_d <- median(diff_wp)
wp_low_d <- quantile(diff_wp, 0.025)
wp_up_d <- quantile(diff_wp, 0.975)
wp0.00 <- mean(diff_wp > 0)
wp0.02 <- mean(diff_wp > 2)
wp0.05 <- mean(diff_wp > 5)
wp0.10 <- mean(diff_wp > 10)
wp0.15 <- mean(diff_wp > 15)

wp_med_r <- median(ratio_wp)
wp_low_r <- quantile(ratio_wp, 0.025)
wp_up_r <- quantile(ratio_wp, 0.975)
wp1.00 <- mean(ratio_wp > 1)
wp1.10 <- mean(ratio_wp > 1.1)
wp1.25 <- mean(ratio_wp > 1.25)
wp1.50 <- mean(ratio_wp > 1.5)
wp2.00 <- mean(ratio_wp > 2)


#######################################################################
# Make benefit tables

# Table on risk difference scale
df <- data.frame("Prior Belief" = c("Optimistic", "", "",
                                    "Neutral", "", "",
                                    "Pessimistic", "", ""),
                 "Strength" = rep(c("Weak", "Moderate", "Strong"), 3),
                 "Median" = c(paste(round(wo_med_d, 1), " (",
                                    round(wo_low_d, 1), ", ",
                                    round(wo_up_d, 1), ")", sep = ""),
                              paste(round(mo_med_d, 1), " (",
                                    round(mo_low_d, 1), ", ",
                                    round(mo_up_d, 1), ")", sep = ""),
                              paste(round(so_med_d, 1), " (",
                                    round(so_low_d, 1), ", ",
                                    round(so_up_d, 1), ")", sep = ""),
                              paste(round(wn_med_d, 1), " (",
                                    round(wn_low_d, 1), ", ",
                                    round(wn_up_d, 1), ")", sep = ""),
                              paste(round(mn_med_d, 1), " (",
                                    round(mn_low_d, 1), ", ",
                                    round(mn_up_d, 1), ")", sep = ""),
                              paste(round(sn_med_d, 1), " (",
                                    round(sn_low_d, 1), ", ",
                                    round(sn_up_d, 1), ")", sep = ""),
                              paste(round(wp_med_d, 1), " (",
                                    round(wp_low_d, 1), ", ",
                                    round(wp_up_d, 1), ")", sep = ""),
                              paste(round(mp_med_d, 1), " (",
                                    round(mp_low_d, 1), ", ",
                                    round(mp_up_d, 1), ")", sep = ""),
                              paste(round(sp_med_d, 1), " (",
                                    round(sp_low_d, 1), ", ",
                                    round(sp_up_d, 1), ")", sep = "")),
                 ">0" = 100*round(c(wo0.00, mo0.00, so0.00,
                                    wn0.00, mn0.00, sn0.00,
                                    wp0.00, mp0.00, sp0.00), 2),
                 ">2" = 100*round(c(wo0.02, mo0.02, so0.02,
                                    wn0.02, mn0.02, sn0.02,
                                    wp0.02, mp0.02, sp0.02), 2),
                 ">5" = 100*round(c(wo0.05, mo0.05, so0.05,
                                    wn0.05, mn0.05, sn0.05,
                                    wp0.05, mp0.05, sp0.05), 2),
                 ">10" = 100*round(c(wo0.10, mo0.10, so0.10,
                                     wn0.10, mn0.10, sn0.10,
                                     wp0.10, mp0.10, sp0.10), 2),
                 ">15" = 100*round(c(wo0.15, mo0.15, so0.15,
                                     wn0.15, mn0.15, sn0.15,
                                     wp0.15, mp0.15, sp0.15), 2))

ft <- flextable(df) %>%
  theme_zebra()
ft <- width(ft, j = 3, width = 1.5)

# Export table to Word
read_docx() %>% 
  body_add_par("Survival Benefit Table (RD Scale)") %>% 
  body_add_flextable(value = ft) %>% 
  print(target = "survival_benefit_table_RD.docx")


# Now do table on relative risk scale
df2 <- data.frame("Prior Belief" = c("Optimistic", "", "",
                                     "Neutral", "", "",
                                     "Pessimistic", "", ""),
                  "Strength" = rep(c("Weak", "Moderate", "Strong"), 3),
                  "Median" = c(paste(round(wo_med_r, 2), " (",
                                     round(wo_low_r, 2), ", ",
                                     round(wo_up_r, 2), ")", sep = ""),
                               paste(round(mo_med_r, 2), " (",
                                     round(mo_low_r, 2), ", ",
                                     round(mo_up_r, 2), ")", sep = ""),
                               paste(round(so_med_r, 2), " (",
                                     round(so_low_r, 2), ", ",
                                     round(so_up_r, 2), ")", sep = ""),
                               paste(round(wn_med_r, 2), " (",
                                     round(wn_low_r, 2), ", ",
                                     round(wn_up_r, 2), ")", sep = ""),
                               paste(round(mn_med_r, 2), " (",
                                     round(mn_low_r, 2), ", ",
                                     round(mn_up_r, 2), ")", sep = ""),
                               paste(round(sn_med_r, 2), " (",
                                     round(sn_low_r, 2), ", ",
                                     round(sn_up_r, 2), ")", sep = ""),
                               paste(round(wp_med_r, 2), " (",
                                     round(wp_low_r, 2), ", ",
                                     round(wp_up_r, 2), ")", sep = ""),
                               paste(round(mp_med_r, 2), " (",
                                     round(mp_low_r, 2), ", ",
                                     round(mp_up_r, 2), ")", sep = ""),
                               paste(round(sp_med_r, 2), " (",
                                     round(sp_low_r, 2), ", ",
                                     round(sp_up_r, 2), ")", sep = "")),
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

# Export table to Word (man. pull flat prior values from flat-prior-figure.R)
read_docx() %>% 
  body_add_par("Benefit Table (RR Scale)") %>% 
  body_add_flextable(value = ft2) %>% 
  print(target = "benefit_table_RR.docx")