# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(flextable)
library(insight)
library(bayestestR)
library(bayesplot)
library(cowplot)
#library(brms)
library(broom)
library(parameters)
library(sanon)
library(rstanarm)

# Load trial data
dat <- read.csv("/Users/blette/Downloads/outcomes.csv")
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
dat_primary$AgeFactor2 <- 1*(dat_primary$AgeGroup == 2)
dat_primary$AgeFactor3 <- 1*(dat_primary$AgeGroup == 3)

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


#############################################################################
# Run analysis for evidence-based priors

# Prior using Grandfelt 2021 meta-analysis and 50% weighting
grandfelt_mean <- log(1.2)
grandfelt_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.85),
                           priormean = grandfelt_mean)
dw <- 0.5
gf_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeFactor2 + AgeFactor3,
                   data = dat_primary, prior_intercept = normal(0, 100),
                   prior = normal(location = c(grandfelt_mean, 0, 0),
                                  scale = c(sqrt(grandfelt_sd^2 / dw),
                                            100, 100)),
                   chains = 8, family = binomial(link = "log"),
                   seed = 1234, iter = 10000)

pred1_gf <- posterior_epred(gf_mod,
                            newdata = mutate(gf_mod$data, Trt = 1))
pred0_gf <- posterior_epred(gf_mod,
                            newdata = mutate(gf_mod$data, Trt = 0))

ratio_gf <- apply(pred1_gf, 1, mean) / apply(pred0_gf, 1, mean)
diff_gf <- 100*(apply(pred1_gf, 1, mean) - apply(pred0_gf, 1, mean))

# Prior using TTM and TTM2 and 50% weighting
ttm_mean <- log(1)
ttm_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.9),
                     priormean = ttm_mean)
dw <- 0.5
ttm_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeFactor2 + AgeFactor3,
                    data = dat_primary, prior_intercept = normal(0, 100),
                    prior = normal(location = c(ttm_mean, 0, 0),
                                   scale = c(sqrt(ttm_sd^2 / dw), 100, 100)),
                    chains = 8, family = binomial(link = "log"),
                    seed = 1234, iter = 10000)

pred1_ttm <- posterior_epred(ttm_mod,
                             newdata = mutate(ttm_mod$data, Trt = 1))
pred0_ttm <- posterior_epred(ttm_mod,
                             newdata = mutate(ttm_mod$data, Trt = 0))

ratio_ttm <- apply(pred1_ttm, 1, mean) / apply(pred0_ttm, 1, mean)
diff_ttm <- 100*(apply(pred1_ttm, 1, mean) - apply(pred0_ttm, 1, mean))

# Prior using Hyperion with 50% weighting
hyp_mean <- log(1.8)
hyp_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.7),
                     priormean = hyp_mean)
dw <- 0.5
hyp_mod <- stan_glm(PrimaryEndpoint ~ Trt + AgeFactor2 + AgeFactor3,
                    data = dat_primary, prior_intercept = normal(0, 100),
                    prior = normal(location = c(hyp_mean, 0, 0),
                                   scale = c(sqrt(hyp_sd^2 / dw), 100, 100)),
                    chains = 8, family = binomial(link = "log"),
                    seed = 1234, iter = 10000)

pred1_hyp <- posterior_epred(hyp_mod,
                             newdata = mutate(hyp_mod$data, Trt = 1))
pred0_hyp <- posterior_epred(hyp_mod,
                             newdata = mutate(hyp_mod$data, Trt = 0))

ratio_hyp <- apply(pred1_hyp, 1, mean) / apply(pred0_hyp, 1, mean)
diff_hyp <- 100*(apply(pred1_hyp, 1, mean) - apply(pred0_hyp, 1, mean))


#######################################################################
# Prepare to make RD figure
# Helper function to get approximate prior on RD scale from log(OR) scale
approxPrior <- function(num, baserisk, mean, sd) {
  logORprior <- rnorm(num, mean, sd)
  intermed <- exp(logORprior)*baserisk / (1 - baserisk)
  return(100*(intermed / (1 + intermed) - baserisk))
}

# Make data set with densities for each of the 3 priors
temp_plot <- ggplot(data.frame(diff_sn), aes(x = diff_sn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_sn),
                                      apply(pred0_sn, 1, mean), 0,
                                      sn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn2 <- p$data[[1]][, c(1, 2)]
plot_data_sn2$Distribution <- "Prior"

plot_data_sn <- rbind(plot_data_sn, plot_data_sn2)
plot_data_sn$barriers <- 0
plot_data_sn$barriers[plot_data_sn$x < 0] <- 1
plot_data_sn$strength <- "Strong"
plot_data_sn$belief <- "Neutral"
plot_data_sn$label <- paste0("ABD = ",
                             round(median(diff_sn), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_sn, 0.025), 1), ", ",
                             round(quantile(diff_sn, 0.975), 1), ")")