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
# Helper function to get approximate prior on RD scale from log(RR) scale
approxPrior <- function(num, baserisk, mean, sd) {
  logRRprior <- rnorm(num, mean, sd)
  intermed <- exp(logRRprior)*baserisk
  return(100*(intermed - baserisk))
}

# Make data set with densities for each of the 3 priors
temp_plot <- ggplot(data.frame(diff_gf), aes(x = diff_gf)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_gf <- p$data[[1]][, c(1, 2)]
plot_data_gf$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_gf),
                                      apply(pred0_gf, 1, mean), grandfelt_mean,
                                      sqrt(grandfelt_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_gf2 <- p$data[[1]][, c(1, 2)]
plot_data_gf2$Distribution <- "Prior"

plot_data_gf <- rbind(plot_data_gf, plot_data_gf2)
plot_data_gf$barriers <- 0
plot_data_gf$barriers[plot_data_gf$x < 0] <- 1
plot_data_gf$EBP <- "Grandfelt"
plot_data_gf$label <- paste0("ABD = ",
                             round(median(diff_gf), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_gf, 0.025), 1), ", ",
                             round(quantile(diff_gf, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_ttm), aes(x = diff_ttm)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_ttm <- p$data[[1]][, c(1, 2)]
plot_data_ttm$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_ttm),
                                      apply(pred0_ttm, 1, mean), ttm_mean,
                                      sqrt(ttm_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_ttm2 <- p$data[[1]][, c(1, 2)]
plot_data_ttm2$Distribution <- "Prior"

plot_data_ttm <- rbind(plot_data_ttm, plot_data_ttm2)
plot_data_ttm$barriers <- 0
plot_data_ttm$barriers[plot_data_ttm$x < 0] <- 1
plot_data_ttm$EBP <- "TTM"
plot_data_ttm$label <- paste0("ABD = ",
                              round(median(diff_ttm), 1),
                              "\n", "95% CI: (",
                              round(quantile(diff_ttm, 0.025), 1), ", ",
                              round(quantile(diff_ttm, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_hyp), aes(x = diff_hyp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_hyp <- p$data[[1]][, c(1, 2)]
plot_data_hyp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_hyp),
                                      apply(pred0_hyp, 1, mean), hyp_mean,
                                      sqrt(hyp_sd^2 / dw))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_hyp2 <- p$data[[1]][, c(1, 2)]
plot_data_hyp2$Distribution <- "Prior"

plot_data_hyp <- rbind(plot_data_hyp, plot_data_hyp2)
plot_data_hyp$barriers <- 0
plot_data_hyp$barriers[plot_data_hyp$x < 0] <- 1
plot_data_hyp$EBP <- "Hyperion"
plot_data_hyp$label <- paste0("ABD = ",
                              round(median(diff_hyp), 1),
                              "\n", "95% CI: (",
                              round(quantile(diff_hyp, 0.025), 1), ", ",
                              round(quantile(diff_hyp, 0.975), 1), ")")

plot_data <- rbind(plot_data_gf, plot_data_ttm, plot_data_hyp)

f_labels <-
  data.frame(EBP = c("Grandfelt", "TTM", "Hyperion"),
             Distribution = rep("Posterior", 3),
             barriers = rep(0, 9),
             label = c(paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_gf)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_gf,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_gf,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_ttm)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_ttm,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_ttm,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_hyp)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_hyp,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_hyp,
                                                       0.975)), ")")))

# Make figure for RD scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_wrap(~EBP, scales = "free") +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(labels = seq(-20, 20, by = 10),
                     breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  coord_cartesian(xlim = c(-25, 30),
                  ylim = c(0, 0.4)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_bw() +
  scale_fill_brewer(guide = "none") +
  geom_vline(xintercept = 0, color = "black",
             linetype = 1) +
  geom_line() +
  geom_text(data = f_labels, size = 2,
            aes(x = 20, y = 0.16, label = label))

# Output figure
pdf("ebp-neuro-posterior-density-RD-figure.pdf", width = 8.5, height = 6)
fig1
dev.off()

#######################################################################
# Now make figure on relative risk scale

# Plot densities for each of the 3 priors
temp_plot <- ggplot(data.frame(ratio_gf), aes(x = ratio_gf)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_gf <- p$data[[1]][, c(1, 2)]
plot_data_gf$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = exp(rnorm(length(ratio_gf), grandfelt_mean,
                                    sqrt(grandfelt_sd^2 / dw)))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_gf2 <- p$data[[1]][, c(1, 2)]
plot_data_gf2$Distribution <- "Prior"

plot_data_gf <- rbind(plot_data_gf, plot_data_gf2)
plot_data_gf$barriers <- 0
plot_data_gf$barriers[plot_data_gf$x < 1] <- 1
plot_data_gf$EBP <- "Grandfelt"
plot_data_gf$label <- paste0("RB = ",
                             round(median(ratio_gf), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_gf, 0.025), 2), ", ",
                             round(quantile(ratio_gf, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_ttm), aes(x = ratio_ttm)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_ttm <- p$data[[1]][, c(1, 2)]
plot_data_ttm$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = exp(rnorm(length(ratio_ttm), ttm_mean,
                                    sqrt(ttm_sd^2 / dw)))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_ttm2 <- p$data[[1]][, c(1, 2)]
plot_data_ttm2$Distribution <- "Prior"

plot_data_ttm <- rbind(plot_data_ttm, plot_data_ttm2)
plot_data_ttm$barriers <- 0
plot_data_ttm$barriers[plot_data_ttm$x < 1] <- 1
plot_data_ttm$EBP <- "TTM"
plot_data_ttm$label <- paste0("RB = ",
                              round(median(ratio_ttm), 2),
                              "\n", "95% CI: (",
                              round(quantile(ratio_ttm, 0.025), 2), ", ",
                              round(quantile(ratio_ttm, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_hyp), aes(x = ratio_hyp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_hyp <- p$data[[1]][, c(1, 2)]
plot_data_hyp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = exp(rnorm(length(ratio_hyp), hyp_mean,
                                    sqrt(hyp_sd^2 / dw)))),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_hyp2 <- p$data[[1]][, c(1, 2)]
plot_data_hyp2$Distribution <- "Prior"

plot_data_hyp <- rbind(plot_data_hyp, plot_data_hyp2)
plot_data_hyp$barriers <- 0
plot_data_hyp$barriers[plot_data_hyp$x < 1] <- 1
plot_data_hyp$EBP <- "Hyperion"
plot_data_hyp$label <- paste0("RB = ",
                              round(median(ratio_hyp), 2),
                              "\n", "95% CI: (",
                              round(quantile(ratio_hyp, 0.025), 2), ", ",
                              round(quantile(ratio_hyp, 0.975), 2), ")")

plot_data <- rbind(plot_data_gf, plot_data_ttm, plot_data_hyp)

f_labels <-
  data.frame(EBP = c("Grandfelt", "TTM", "Hyperion"),
             Distribution = rep("Posterior", 3),
             barriers = rep(0, 9),
             label = c(paste0("Median RB = ",
                              sprintf("%.1f", median(ratio_gf)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(ratio_gf,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(ratio_gf,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.1f", median(ratio_ttm)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(ratio_ttm,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(ratio_ttm,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.1f", median(ratio_hyp)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(ratio_hyp,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(ratio_hyp,
                                                       0.975)), ")")))

# Make figure for RR scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_wrap(~EBP, scales = "free") +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(trans = "log", labels = seq(0.5, 3, by = 0.5),
                     breaks = seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  coord_cartesian(xlim = c(0.5, 3.4),
                  ylim = c(0, 5.75)) +
  geom_ribbon(aes(ymin=0, ymax=y, fill=factor(barriers)),
              show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "black",
             linetype = 1) +
  theme_bw() +
  scale_fill_brewer(guide = "none") +
  geom_vline(xintercept = 1, color = "black",
             linetype = 1) +
  geom_line() +
  geom_text(data = f_labels, size = 2,
            aes(x = 2.25, y = 3, label = label))

# Output figure
pdf("ebp-neuro-posterior-density-RR-figure.pdf", width = 8.5, height = 6)
fig2
dev.off()


##########################################################################
# Run model diagnostics

# Borrow diagnostics function from COVID Steroid 2 Trial analysis
# MCMC diagnostics
diag_fit <- function(fit, extra_groups = NULL, bars = TRUE) {
  
  # Print fit (includes Rhats and bulk and tail ESS)
  print(summary(fit))
  if (any(rhat(fit) > 1.01)) {
    warning("One or more Rhats > 1.01")
  } else {
    message("All Rhats <= 1.01")
  }
  
  #nam <- names(fit$fit)
  #walk(nam, ~{print(mcmc_trace(fit, pars = .x))})
  print(mcmc_trace(fit))
  #walk(nam, ~{print(mcmc_dens_overlay(fit, pars = .x))})
  print(mcmc_dens_overlay(fit))
  
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
  
  print(loo(fit, cores = 4)) # RELOO IF NECESSARY
  invisible(fit)
  
}

diag_fit(gf_mod)
diag_fit(ttm_mod)
diag_fit(hyp_mod)
