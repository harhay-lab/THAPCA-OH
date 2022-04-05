# Clear data and load libraries
rm(list = ls())
library(tidyverse)
library(flextable)
library(insight)
library(bayestestR)
library(bayesplot)
library(cowplot)
library(brms)
library(broom)
library(parameters)
library(sanon)
#library(rstanarm)

# Load trial data
dat <- read.csv("/Users/blette/Downloads/outcomes.csv")
dat_primary <- dat[dat$Primary == 1 & dat$PrimaryEndpoint != 97, ]
dat_primary$Trt <- abs(2 - dat_primary$TreatRand)
dat_primary$AgeFactor <- as.factor(dat_primary$AgeGroup)

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
# Analysis

# Optimistic priors
# Strong optimistic prior
o_mean <- log(1.25)
so_sd <- getPriorSD(propbelow = 0.05, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(so_sd)
so_prior <- prior(normal(o_mean, so_sd), class = "b", coef = Trt)
so_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = so_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)
  
pred1_so <- posterior_epred(so_mod,
                            newdata = mutate(so_mod$data, Trt = 1))
pred0_so <- posterior_epred(so_mod,
                            newdata = mutate(so_mod$data, Trt = 0))

ratio_so <- apply(pred1_so, 1, mean) / apply(pred0_so, 1, mean)
diff_so <- 100*(apply(pred1_so, 1, mean) - apply(pred0_so, 1, mean))

# Moderate optimistic prior
mo_sd <- getPriorSD(propbelow = 0.15, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(mo_sd)
mo_prior <- prior(normal(o_mean, mo_sd), class = "b", coef = Trt)
mo_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = mo_prior, stanvars = stanvars, iter = 10000,
              family = "bernoulli", seed = 1234, chains = 8)

pred1_mo <- posterior_epred(mo_mod,
                            newdata = mutate(mo_mod$data,Trt = 1))
pred0_mo <- posterior_epred(mo_mod,
                            newdata = mutate(mo_mod$data, Trt = 0))

ratio_mo <- apply(pred1_mo, 1, mean) / apply(pred0_mo, 1, mean)
diff_mo <- 100*(apply(pred1_mo, 1, mean) - apply(pred0_mo, 1, mean))

# Weak optimistic prior
wo_sd <- getPriorSD(propbelow = 0.3, belowcutoff = 0,
                    priormean = o_mean)
stanvars <- stanvar(o_mean) + stanvar(wo_sd)
wo_prior <- prior(normal(o_mean, wo_sd), class = "b", coef = Trt)
wo_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = wo_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wo <- posterior_epred(wo_mod,
                            newdata = mutate(wo_mod$data, Trt = 1))
pred0_wo <- posterior_epred(wo_mod,
                            newdata = mutate(wo_mod$data, Trt = 0))

ratio_wo <- apply(pred1_wo, 1, mean) / apply(pred0_wo, 1, mean)
diff_wo <- 100*(apply(pred1_wo, 1, mean) - apply(pred0_wo, 1, mean))

# Neutral priors
# Strong neutral prior
sn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(1/1.5),
                    priormean = 0)
stanvars <- stanvar(sn_sd)
sn_prior <- prior(normal(0, sn_sd), class = "b", coef = Trt)
sn_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = sn_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_sn <- posterior_epred(sn_mod,
                            newdata = mutate(sn_mod$data, Trt = 1))
pred0_sn <- posterior_epred(sn_mod,
                            newdata = mutate(sn_mod$data, Trt = 0))

ratio_sn <- apply(pred1_sn, 1, mean) / apply(pred0_sn, 1, mean)
diff_sn <- 100*(apply(pred1_sn, 1, mean) - apply(pred0_sn, 1, mean))

# Moderate neutral prior
mn_sd <- getPriorSD(propbelow = 0.025, belowcutoff = log(0.5),
                    priormean = 0)
stanvars <- stanvar(mn_sd)
mn_prior <- prior(normal(0, mn_sd), class = "b", coef = Trt)
mn_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = mn_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_mn <- posterior_epred(mn_mod,
                            newdata = mutate(mn_mod$data, Trt = 1))
pred0_mn <- posterior_epred(mn_mod,
                            newdata = mutate(mn_mod$data, Trt = 0))

ratio_mn <- apply(pred1_mn, 1, mean) / apply(pred0_mn, 1, mean)
diff_mn <- 100*(apply(pred1_mn, 1, mean) - apply(pred0_mn, 1, mean))

# Weak neutral prior
wn_sd <- 3
wn_prior <- prior(normal(0, 3), class = "b", coef = Trt)
wn_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = wn_prior, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wn <- posterior_epred(wn_mod,
                            newdata = mutate(wn_mod$data, Trt = 1))
pred0_wn <- posterior_epred(wn_mod,
                            newdata = mutate(wn_mod$data, Trt = 0))

ratio_wn <- apply(pred1_wn, 1, mean) / apply(pred0_wn, 1, mean)
diff_wn <- 100*(apply(pred1_wn, 1, mean) - apply(pred0_wn, 1, mean))

# Pessimistic priors
# Strong pessimistic prior
p_mean <- log(1/1.25)
sp_sd <- so_sd
stanvars <- stanvar(p_mean) + stanvar(sp_sd)
sp_prior <- prior(normal(p_mean, sp_sd), class = "b", coef = Trt)
sp_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = sp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_sp <- posterior_epred(sp_mod,
                            newdata = mutate(sp_mod$data, Trt = 1))
pred0_sp <- posterior_epred(sp_mod,
                            newdata = mutate(sp_mod$data, Trt = 0))

ratio_sp <- apply(pred1_sp, 1, mean) / apply(pred0_sp, 1, mean)
diff_sp <- 100*(apply(pred1_sp, 1, mean) - apply(pred0_sp, 1, mean))

# Moderate pessimistic prior
mp_sd <- mo_sd
stanvars <- stanvar(p_mean) + stanvar(mp_sd)
mp_prior <- prior(normal(p_mean, mp_sd), class = "b", coef = Trt)
mp_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = mp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_mp <- posterior_epred(mp_mod,
                            newdata = mutate(mp_mod$data, Trt = 1))
pred0_mp <- posterior_epred(mp_mod,
                            newdata = mutate(mp_mod$data, Trt = 0))

ratio_mp <- apply(pred1_mp, 1, mean) / apply(pred0_mp, 1, mean)
diff_mp <- 100*(apply(pred1_mp, 1, mean) - apply(pred0_mp, 1, mean))

# Weak pessimistic prior
wp_sd <- wo_sd
stanvars <- stanvar(p_mean) + stanvar(wp_sd)
wp_prior <- prior(normal(p_mean, wp_sd), class = "b", coef = Trt)
wp_mod <- brm(PrimaryEndpoint ~ Trt + AgeFactor, data = dat_primary,
              prior = wp_prior, stanvars = stanvars, chains = 8,
              family = "bernoulli", seed = 1234, iter = 10000)

pred1_wp <- posterior_epred(wp_mod,
                            newdata = mutate(wp_mod$data, Trt = 1))
pred0_wp <- posterior_epred(wp_mod,
                            newdata = mutate(wp_mod$data, Trt = 0))

ratio_wp <- apply(pred1_wp, 1, mean) / apply(pred0_wp, 1, mean)
diff_wp <- 100*(apply(pred1_wp, 1, mean) - apply(pred0_wp, 1, mean))


#######################################################################
# Prepare to make RD figure
# Helper function to get approximate prior on RD scale from log(OR) scale
approxPrior <- function(num, baserisk, mean, sd) {
  logORprior <- rnorm(num, mean, sd)
  intermed <- exp(logORprior)*baserisk / (1 - baserisk)
  return(100*(intermed / (1 + intermed) - baserisk))
}

# Can I pull baseline risk for each model output line?
#hist(approxPrior(length(diff_wp), apply(pred0, 1, mean), p_mean, wp_sd))
#hist(approxPrior(length(diff_wn), apply(pred0, 1, mean), 0, 3))

# Make data set with densities for each of the 9 priors
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

temp_plot <- ggplot(data.frame(diff_mn), aes(x = diff_mn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn <- p$data[[1]][, c(1, 2)]
plot_data_mn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_mn),
                                      apply(pred0_mn, 1, mean), 0,
                                      mn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn2 <- p$data[[1]][, c(1, 2)]
plot_data_mn2$Distribution <- "Prior"

plot_data_mn <- rbind(plot_data_mn, plot_data_mn2)
plot_data_mn$barriers <- 0
plot_data_mn$barriers[plot_data_mn$x < 0] <- 1
plot_data_mn$strength <- "Moderate"
plot_data_mn$belief <- "Neutral"
plot_data_mn$label <- paste0("ABD = ",
                             round(median(diff_mn), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_mn, 0.025), 1), ", ",
                             round(quantile(diff_mn, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_wn), aes(x = diff_wn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn <- p$data[[1]][, c(1, 2)]
plot_data_wn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_wn),
                                      apply(pred0_wn, 1, mean), 0,
                                      wn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn2 <- p$data[[1]][, c(1, 2)]
plot_data_wn2$Distribution <- "Prior"

plot_data_wn <- rbind(plot_data_wn, plot_data_wn2)
plot_data_wn$barriers <- 0
plot_data_wn$barriers[plot_data_wn$x < 0] <- 1
plot_data_wn$strength <- " Weak"
plot_data_wn$belief <- "Neutral"
plot_data_wn$label <- paste0("ABD = ",
                             round(median(diff_wn), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_wn, 0.025), 1), ", ",
                             round(quantile(diff_wn, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_so), aes(x = diff_so)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so <- p$data[[1]][, c(1, 2)]
plot_data_so$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_so),
                                      apply(pred0_so, 1, mean), o_mean,
                                      so_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so2 <- p$data[[1]][, c(1, 2)]
plot_data_so2$Distribution <- "Prior"

plot_data_so <- rbind(plot_data_so, plot_data_so2)
plot_data_so$barriers <- 0
plot_data_so$barriers[plot_data_so$x < 0] <- 1
plot_data_so$strength <- "Strong"
plot_data_so$belief <- " Optimistic"
plot_data_so$label <- paste0("ABD = ",
                             round(median(diff_so), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_so, 0.025), 1), ", ",
                             round(quantile(diff_so, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_mo), aes(x = diff_mo)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo <- p$data[[1]][, c(1, 2)]
plot_data_mo$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_mo),
                                      apply(pred0_mo, 1, mean), o_mean,
                                      mo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo2 <- p$data[[1]][, c(1, 2)]
plot_data_mo2$Distribution <- "Prior"

plot_data_mo <- rbind(plot_data_mo, plot_data_mo2)
plot_data_mo$barriers <- 0
plot_data_mo$barriers[plot_data_mo$x < 0] <- 1
plot_data_mo$strength <- "Moderate"
plot_data_mo$belief <- " Optimistic"
plot_data_mo$label <- paste0("ABD = ",
                             round(median(diff_mo), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_mo, 0.025), 1), ", ",
                             round(quantile(diff_mo, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_wo), aes(x = diff_wo)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo <- p$data[[1]][, c(1, 2)]
plot_data_wo$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_wo),
                                      apply(pred0_wo, 1, mean), o_mean,
                                      wo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo2 <- p$data[[1]][, c(1, 2)]
plot_data_wo2$Distribution <- "Prior"

plot_data_wo <- rbind(plot_data_wo, plot_data_wo2)
plot_data_wo$barriers <- 0
plot_data_wo$barriers[plot_data_wo$x < 0] <- 1
plot_data_wo$strength <- " Weak"
plot_data_wo$belief <- " Optimistic"
plot_data_wo$label <- paste0("ABD = ",
                             round(median(diff_wo), 1),
                             "\n", "95% CrI: (",
                             round(quantile(diff_wo, 0.025), 1), ", ",
                             round(quantile(diff_wo, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_sp), aes(x = diff_sp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp <- p$data[[1]][, c(1, 2)]
plot_data_sp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_sp),
                                      apply(pred0_sp, 1, mean), p_mean,
                                      sp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp2 <- p$data[[1]][, c(1, 2)]
plot_data_sp2$Distribution <- "Prior"

plot_data_sp <- rbind(plot_data_sp, plot_data_sp2)
plot_data_sp$barriers <- 0
plot_data_sp$barriers[plot_data_sp$x < 0] <- 1
plot_data_sp$strength <- "Strong"
plot_data_sp$belief <- "Pessimistic"
plot_data_sp$label <- paste0("ABD = ",
                             round(median(diff_sp), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_sp, 0.025), 1), ", ",
                             round(quantile(diff_sp, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_mp), aes(x = diff_mp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp <- p$data[[1]][, c(1, 2)]
plot_data_mp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_mp),
                                      apply(pred0_mp, 1, mean), p_mean,
                                      mp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp2 <- p$data[[1]][, c(1, 2)]
plot_data_mp2$Distribution <- "Prior"

plot_data_mp <- rbind(plot_data_mp, plot_data_mp2)
plot_data_mp$barriers <- 0
plot_data_mp$barriers[plot_data_mp$x < 0] <- 1
plot_data_mp$strength <- "Moderate"
plot_data_mp$belief <- "Pessimistic"
plot_data_mp$label <- paste0("ABD = ",
                             round(median(diff_mp), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_mp, 0.025), 1), ", ",
                             round(quantile(diff_mp, 0.975), 1), ")")

temp_plot <- ggplot(data.frame(diff_wp), aes(x = diff_wp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp <- p$data[[1]][, c(1, 2)]
plot_data_wp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior(length(diff_wp),
                                      apply(pred0_wp, 1, mean), p_mean,
                                      wp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp2 <- p$data[[1]][, c(1, 2)]
plot_data_wp2$Distribution <- "Prior"

plot_data_wp <- rbind(plot_data_wp, plot_data_wp2)
plot_data_wp$barriers <- 0
plot_data_wp$barriers[plot_data_wp$x < 0] <- 1
plot_data_wp$strength <- " Weak"
plot_data_wp$belief <- "Pessimistic"
plot_data_wp$label <- paste0("ABD = ",
                             round(median(diff_wp), 1),
                             "\n", "95% CI: (",
                             round(quantile(diff_wp, 0.025), 1), ", ",
                             round(quantile(diff_wp, 0.975), 1), ")")

plot_data <- rbind(plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sp, plot_data_mp, plot_data_wp)

f_labels <-
  data.frame(belief = c(" Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic"),
             strength = c("Strong", "Moderate", " Weak",
                          "Moderate", " Weak", "Strong",
                          " Weak", "Strong", "Moderate"),
             Distribution = rep("Posterior", 9),
             barriers = rep(0, 9),
             label = c(paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_so)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_so,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_so,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_mn)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_mn,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_mn,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_wp)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_wp,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_wp,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_mo)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_mo,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_mo,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_wn)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_wn,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_wn,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_sp)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_sp,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_sp,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_wo)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_wo,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_wo,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_sn)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_sn,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_sn,
                                                       0.975)), ")"),
                       paste0("Median ABD = ",
                              sprintf("%.1f", median(diff_mp)),
                              "\n", "95% CrI: (",
                              sprintf("%.1f", quantile(diff_mp,
                                                       0.025)), ", ",
                              sprintf("%.1f", quantile(diff_mp,
                                                       0.975)), ")")))

# Make figure for RD scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig1 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_grid(belief ~ strength, scales = "free", switch = "y") +
  labs(x = "Absolute Benefit Difference (%)", y = "Density") +
  scale_x_continuous(labels = seq(-20, 20, by = 10),
                     breaks = seq(-20, 20, by = 10)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  coord_cartesian(xlim = c(-25, 30),
                  ylim = c(0, 0.25)) +
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
pdf("posterior-density-RD-figure.pdf", width = 8.5, height = 6)
fig1
dev.off()

#######################################################################
# Now make figure on relative risk scale

# Helper function to get approximate prior on RR scale from log(OR) scale
approxPrior2 <- function(num, baserisk, mean, sd) {
  logORprior <- rnorm(num, mean, sd)
  intermed <- exp(logORprior)*baserisk / (1 - baserisk)
  return(intermed / (1 + intermed) / baserisk)
}

# Plot densities for each of the 9 priors
temp_plot <- ggplot(data.frame(ratio_sn), aes(x = ratio_sn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn <- p$data[[1]][, c(1, 2)]
plot_data_sn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_sn),
                                       apply(pred0_sn, 1, mean), 0,
                                       sn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sn2 <- p$data[[1]][, c(1, 2)]
plot_data_sn2$Distribution <- "Prior"

plot_data_sn <- rbind(plot_data_sn, plot_data_sn2)
plot_data_sn$barriers <- 0
plot_data_sn$barriers[plot_data_sn$x < 1] <- 1
plot_data_sn$strength <- "Strong"
plot_data_sn$belief <- "Neutral"
plot_data_sn$label <- paste0("RB = ",
                             round(median(ratio_sn), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_sn, 0.025), 2), ", ",
                             round(quantile(ratio_sn, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_mn), aes(x = ratio_mn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn <- p$data[[1]][, c(1, 2)]
plot_data_mn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_mn),
                                       apply(pred0_mn, 1, mean), 0,
                                       mn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mn2 <- p$data[[1]][, c(1, 2)]
plot_data_mn2$Distribution <- "Prior"

plot_data_mn <- rbind(plot_data_mn, plot_data_mn2)
plot_data_mn$barriers <- 0
plot_data_mn$barriers[plot_data_mn$x < 1] <- 1
plot_data_mn$strength <- "Moderate"
plot_data_mn$belief <- "Neutral"
plot_data_mn$label <- paste0("RB = ",
                             round(median(ratio_mn), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_mn, 0.025), 2), ", ",
                             round(quantile(ratio_mn, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_wn), aes(x = ratio_wn)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn <- p$data[[1]][, c(1, 2)]
plot_data_wn$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_wn),
                                       apply(pred0_wn, 1, mean), 0,
                                       wn_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wn2 <- p$data[[1]][, c(1, 2)]
plot_data_wn2$Distribution <- "Prior"

plot_data_wn <- rbind(plot_data_wn, plot_data_wn2)
plot_data_wn$barriers <- 0
plot_data_wn$barriers[plot_data_wn$x < 1] <- 1
plot_data_wn$strength <- " Weak"
plot_data_wn$belief <- "Neutral"
plot_data_wn$label <- paste0("RB = ",
                             round(median(ratio_wn), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_wn, 0.025), 2), ", ",
                             round(quantile(ratio_wn, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_so), aes(x = ratio_so)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so <- p$data[[1]][, c(1, 2)]
plot_data_so$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_so),
                                       apply(pred0_so, 1, mean), o_mean,
                                       so_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_so2 <- p$data[[1]][, c(1, 2)]
plot_data_so2$Distribution <- "Prior"

plot_data_so <- rbind(plot_data_so, plot_data_so2)
plot_data_so$barriers <- 0
plot_data_so$barriers[plot_data_so$x < 1] <- 1
plot_data_so$strength <- "Strong"
plot_data_so$belief <- " Optimistic"
plot_data_so$label <- paste0("RB = ",
                             round(median(ratio_so), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_so, 0.025), 2), ", ",
                             round(quantile(ratio_so, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_mo), aes(x = ratio_mo)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo <- p$data[[1]][, c(1, 2)]
plot_data_mo$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_mo),
                                       apply(pred0_mo, 1, mean), o_mean,
                                       mo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mo2 <- p$data[[1]][, c(1, 2)]
plot_data_mo2$Distribution <- "Prior"

plot_data_mo <- rbind(plot_data_mo, plot_data_mo2)
plot_data_mo$barriers <- 0
plot_data_mo$barriers[plot_data_mo$x < 1] <- 1
plot_data_mo$strength <- "Moderate"
plot_data_mo$belief <- " Optimistic"
plot_data_mo$label <- paste0("RB = ",
                             round(median(ratio_mo), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_mo, 0.025), 2), ", ",
                             round(quantile(ratio_mo, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_wo), aes(x = ratio_wo)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo <- p$data[[1]][, c(1, 2)]
plot_data_wo$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_wo),
                                       apply(pred0_wo, 1, mean), o_mean,
                                       wo_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wo2 <- p$data[[1]][, c(1, 2)]
plot_data_wo2$Distribution <- "Prior"

plot_data_wo <- rbind(plot_data_wo, plot_data_wo2)
plot_data_wo$barriers <- 0
plot_data_wo$barriers[plot_data_wo$x < 1] <- 1
plot_data_wo$strength <- " Weak"
plot_data_wo$belief <- " Optimistic"
plot_data_wo$label <- paste0("RB = ",
                             round(median(ratio_wo), 2),
                             "\n", "95% CrI: (",
                             round(quantile(ratio_wo, 0.025), 2), ", ",
                             round(quantile(ratio_wo, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_sp), aes(x = ratio_sp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp <- p$data[[1]][, c(1, 2)]
plot_data_sp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_sp),
                                       apply(pred0_sp, 1, mean), p_mean,
                                       sp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_sp2 <- p$data[[1]][, c(1, 2)]
plot_data_sp2$Distribution <- "Prior"

plot_data_sp <- rbind(plot_data_sp, plot_data_sp2)
plot_data_sp$barriers <- 0
plot_data_sp$barriers[plot_data_sp$x < 1] <- 1
plot_data_sp$strength <- "Strong"
plot_data_sp$belief <- "Pessimistic"
plot_data_sp$label <- paste0("RB = ",
                             round(median(ratio_sp), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_sp, 0.025), 2), ", ",
                             round(quantile(ratio_sp, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_mp), aes(x = ratio_mp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp <- p$data[[1]][, c(1, 2)]
plot_data_mp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_mp),
                                       apply(pred0_mp, 1, mean), p_mean,
                                       mp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_mp2 <- p$data[[1]][, c(1, 2)]
plot_data_mp2$Distribution <- "Prior"

plot_data_mp <- rbind(plot_data_mp, plot_data_mp2)
plot_data_mp$barriers <- 0
plot_data_mp$barriers[plot_data_mp$x < 1] <- 1
plot_data_mp$strength <- "Moderate"
plot_data_mp$belief <- "Pessimistic"
plot_data_mp$label <- paste0("RB = ",
                             round(median(ratio_mp), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_mp, 0.025), 2), ", ",
                             round(quantile(ratio_mp, 0.975), 2), ")")

temp_plot <- ggplot(data.frame(ratio_wp), aes(x = ratio_wp)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp <- p$data[[1]][, c(1, 2)]
plot_data_wp$Distribution <- "Posterior"

temp_plot <-
  ggplot(data.frame(pri = approxPrior2(length(ratio_wp),
                                       apply(pred0_wp, 1, mean), p_mean,
                                       wp_sd)),
         aes(x = pri)) +
  geom_density()
p <- ggplot_build(temp_plot)
plot_data_wp2 <- p$data[[1]][, c(1, 2)]
plot_data_wp2$Distribution <- "Prior"

plot_data_wp <- rbind(plot_data_wp, plot_data_wp2)
plot_data_wp$barriers <- 0
plot_data_wp$barriers[plot_data_wp$x < 1] <- 1
plot_data_wp$strength <- " Weak"
plot_data_wp$belief <- "Pessimistic"
plot_data_wp$label <- paste0("RB = ",
                             round(median(ratio_wp), 2),
                             "\n", "95% CI: (",
                             round(quantile(ratio_wp, 0.025), 2), ", ",
                             round(quantile(ratio_wp, 0.975), 2), ")")

plot_data <- rbind(plot_data_sn, plot_data_mn, plot_data_wn,
                   plot_data_so, plot_data_mo, plot_data_wo,
                   plot_data_sp, plot_data_mp, plot_data_wp)

f_labels <-
  data.frame(belief = c(" Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic",
                        " Optimistic", "Neutral", "Pessimistic"),
             strength = c("Strong", "Moderate", " Weak",
                          "Moderate", " Weak", "Strong",
                          " Weak", "Strong", "Moderate"),
             Distribution = rep("Posterior", 9),
             barriers = rep(0, 9),
             label = c(paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_so)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_so,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_so,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_mn)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_mn,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_mn,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_wp)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_wp,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_wp,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_mo)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_mo,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_mo,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_wn)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_wn,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_wn,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_sp)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_sp,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_sp,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_wo)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_wo,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_wo,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_sn)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_sn,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_sn,
                                                       0.975)), ")"),
                       paste0("Median RB = ",
                              sprintf("%.2f", median(ratio_mp)),
                              "\n", "95% CrI: (",
                              sprintf("%.2f", quantile(ratio_mp,
                                                       0.025)), ", ",
                              sprintf("%.2f", quantile(ratio_mp,
                                                       0.975)), ")")))

# Make figure for RR scale
plot_data$barriers[plot_data$Distribution == "Prior"] <- NA
fig2 <- ggplot(plot_data, aes(x = x, y = y, group = barriers,
                              lty = Distribution)) +
  facet_grid(belief ~ strength, scales = "free", switch = "y") +
  labs(x = "Relative Benefit", y = "Density") +
  scale_x_continuous(trans = "log", labels = seq(0.5, 3, by = 0.5),
                     breaks = seq(0.5, 3, by = 0.5)) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  coord_cartesian(xlim = c(0.5, 3.4),
                  ylim = c(0, 4.75)) +
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
pdf("posterior-density-RR-figure.pdf", width = 8.5, height = 6)
fig2
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

diag_fit(so_mod)
diag_fit(mo_mod)
diag_fit(wo_mod)
diag_fit(sn_mod)
diag_fit(mn_mod)
diag_fit(wn_mod)
diag_fit(sp_mod)
diag_fit(mp_mod)
diag_fit(wp_mod)
